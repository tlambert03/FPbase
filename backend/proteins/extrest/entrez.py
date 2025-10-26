from __future__ import annotations

import contextlib
import datetime
import logging
import os
import re
import time
from typing import TYPE_CHECKING, Literal, TypedDict, cast

from Bio import Entrez, SeqIO
from django.core.cache import cache
from habanero import Crossref

logger = logging.getLogger(__name__)

if TYPE_CHECKING:

    class DoiInfo(TypedDict, total=False):
        doi: str | None
        pmid: str | None
        title: str | None
        journal: str | None
        pages: str | None
        volume: str | None
        issue: str | None
        year: int | None
        date: datetime.date | None
        authors: list[dict] | None


Entrez.email = "talley_lambert@hms.harvard.edu"
Entrez.api_key = os.getenv("NCBI_API_KEY", None)


def _parse_crossref(doidict: dict) -> DoiInfo:
    if "message" in doidict:
        doidict = doidict["message"]
    out = {}
    ct = doidict.get("container-title", [])
    out: DoiInfo = {
        "journal": ct[0] if ct else None,
        "title": doidict.get("title"),
        "pages": doidict.get("page", None),
        "volume": doidict.get("volume", None),
        "issue": doidict.get("issue", None),
        "authors": doidict.get("author", None),
        "doi": doidict.get("DOI", None),
    }
    out["title"] = out["title"][0] if len(out["title"]) else None
    dp = doidict.get("published-online", doidict.get("published-print", doidict.get("issued", {}))).get(
        "date-parts", [None]
    )[0]
    dp.append(1) if dp and len(dp) == 1 else None
    dp.append(1) if dp and len(dp) == 2 else None
    out["date"] = datetime.date(*dp) if dp and all(dp) else None
    out["year"] = dp[0] if dp and dp[0] else None
    return out


def _crossref(doi):
    cr = Crossref(mailto="talley.lambert+fpbase@gmail.org")
    response = cr.works(ids=doi)
    # habanero returns a list if doi is a list of len > 1
    # otherwise a single dict
    if isinstance(doi, list | tuple | set) and len(doi) > 1:
        return {x.pop("doi"): x for x in (_parse_crossref(i) for i in response)}
    else:
        return _parse_crossref(response)


def _doi2pmid(doi: str) -> str | None:
    pubmed_record = Entrez.read(Entrez.esearch(db="pubmed", term=doi))
    try:
        return pubmed_record.get("IdList")[0]
    except Exception:
        return None


def _get_pmid_info(pmid: str) -> DoiInfo | None:
    pubmed_record = Entrez.read(Entrez.esummary(db="pubmed", id=pmid, retmode="xml"))
    if len(pubmed_record):
        pubmed_record = pubmed_record[0]
        date = None
        try:
            date = datetime.datetime.strptime(pubmed_record["PubDate"], "%Y %b %d").date()
        except Exception:
            with contextlib.suppress(Exception):
                date = datetime.datetime.strptime(pubmed_record["EPubDate"], "%Y %b %d").date()
        return {
            "doi": pubmed_record.get("DOI", None),
            "title": pubmed_record.get("Title", ""),
            "journal": pubmed_record.get("Source", ""),
            "pages": pubmed_record.get("Pages", ""),
            "volume": pubmed_record.get("Volume", ""),
            "issue": pubmed_record.get("Issue", ""),
            "year": pubmed_record["PubDate"].split()[0],
            "authors": pubmed_record["AuthorList"],
            "date": date,
        }
    return None


def _merge_info(dict1, dict2, exclude=()):
    """existings values in dict2 will overwrite dict1"""
    for key in dict1.keys():
        if key in dict2 and dict2[key] and key not in exclude:
            dict1[key] = dict2[key]
    return dict1


def doi_lookup(doi: str) -> DoiInfo:
    info = _crossref(doi)
    pmid = _doi2pmid(doi)
    if pmid:
        try:
            pinfo = _get_pmid_info(pmid)
            assert pinfo.pop("doi") == doi
            info = _merge_info(pinfo, info)
            info["pmid"] = pmid
        except AssertionError:
            pass
    # get rid of empty values
    for key, val in list(info.items()):
        if not val:
            del info[key]
    return cast("DoiInfo", info)


def get_cached_gbseqs(gbids, max_age=60 * 60 * 24) -> dict[str, tuple[str, float]]:
    """Get sequences for multiple genbank IDs, cached for max_age seconds.

    Returns dict of gbid -> (sequence, timestamp-fetched)
    {
        'LC381432': (
            'MVSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATIGKLTLKLICTTGKLPVPWPTLVTTL...',
            1761089667.316076
        ),
        'LC085679': (
            'MENVRRKTGIQTEMKTKLHMDGMVNGHSFEIKGEGKGSPYEGVQTMKLKVTKGAPLPFSIDILL...',
            1761089667.316076
        ),
        ...
    }
    """
    gbseqs = cache.get("gbseqs", {})
    now = time.time()
    tofetch = [id for id in gbids if id not in gbseqs or gbseqs[id][1] - now >= max_age]
    gbseqs.update({k: (v, now) for k, v in _fetch_gb_seqs(tofetch).items()})
    cache.set("gbseqs", gbseqs, 60 * 60 * 24)
    return gbseqs


# genbank protein accession
gbprotrx = re.compile(r"^[A-Za-z]{3}\d{5}\.?\d?")
# genbank nucleotide accession
gbnucrx = re.compile(r"^([A-Za-z]{2}\d{6}\.?\d?)|([A-Za-z]{1}\d{5}\.?\d?)")
# refseq accession
refseqrx = re.compile("(NC|AC|NG|NT|NW|NZ|NM|NR|XM|XR|NP|AP|XP|YP|ZP)_[0-9]+")


def _check_accession_type(gbid) -> Literal["protein", "nuccore"] | None:
    if gbprotrx.match(gbid):  # it is a protein accession
        return "protein"
    elif gbnucrx.match(gbid):
        return "nuccore"
    elif gbid.upper().startswith(("WP_", "XP_", "YP_", "NP_", "AP_")):
        return "protein"
    return None


def _fetch_gb_seqs(gbids):
    """Retrieve protein sequence for multiple genbank IDs, (regardless of accession type)"""
    prots = []
    nucs = []
    for id in gbids:
        db = _check_accession_type(id)
        if db == "protein":
            prots.append(id)
        elif db == "nuccore":
            nucs.append(id)
        else:
            prots.append(id)
            logger.error("Could not determine accession type for %s", id)
    records = {}
    if len(nucs):
        with Entrez.efetch(db="nuccore", id=nucs, rettype="fasta", retmode="text") as handle:
            records.update(
                {x.id.split(".")[0]: str(x.translate().seq).strip("*") for x in SeqIO.parse(handle, "fasta")}
            )
    if len(prots):
        with Entrez.efetch(db="protein", id=prots, rettype="fasta", retmode="text") as handle:
            records.update({x.id.split(".")[0]: str(x.seq).strip("*") for x in SeqIO.parse(handle, "fasta")})
    return records


def get_gb_seq(gbid: str) -> str | None:
    """Retrieve protein sequence for genbank ID, (regardless of accession type)"""
    database = _check_accession_type(gbid)
    if not database:
        return None

    if database == "protein":
        with Entrez.efetch(db=database, id=gbid, rettype="fasta", retmode="text") as handle:
            record = SeqIO.read(handle, "fasta")
        if hasattr(record, "seq"):
            return str(record.seq)
    else:
        with Entrez.efetch(db=database, id=gbid, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
        return _parse_gbnuc_record(record).get("seq", None)
    return None


def _parse_gbnuc_record(record):
    """
    returns (possible)
        desc: description
        seq: protein sequence
        pmids: [list of pmids]
        organism: name of organism
        gb_prot: protein accession
        gb_nuc:  nucleotide accession
        product: description of protein product
        gene: addition description of the gene (might have protein name)
    """
    D = {}
    annotations = getattr(record, "annotations", None)
    features = getattr(record, "features", None)
    D["desc"] = getattr(record, "description", None)
    if annotations:
        refs = annotations.get("references")
        if refs:
            D["pmids"] = [getattr(r, "pubmed_id", None) for r in refs if getattr(r, "pubmed_id", False)]
        D["organism"] = annotations.get("organism", None)
        nuc = annotations.get("accessions", [])
        if len(nuc):
            D["gb_nuc"] = nuc[0]
    if features:
        CDSs = [f for f in features if f.type == "CDS"]
        if len(CDSs) == 1:
            cds = CDSs[0]
            quals = getattr(cds, "qualifiers", None)
            if quals:
                trans = quals.get("translation", [])
                if len(trans):
                    D["seq"] = trans[0]
                protid = quals.get("protein_id", [])
                if len(protid):
                    D["gb_prot"] = protid[0].split(".")[0]
                gene = quals.get("gene", [])
                if len(gene):
                    D["gene"] = gene[0]
                product = quals.get("product", [])
                if len(product):
                    D["product"] = product[0]
    return D


def get_organism_info(organism_id: int | str) -> dict | None:
    pubmed_record = Entrez.read(Entrez.esummary(db="taxonomy", id=organism_id, retmode="xml"))
    if not pubmed_record or len(pubmed_record) == 0:
        return None
    r = pubmed_record[0]
    return {
        "scientific_name": r["ScientificName"],
        "division": r["Division"],
        "common_name": r["CommonName"],
        "species": r["Species"],
        "genus": r["Genus"],
        "rank": r["Rank"],
    }

import contextlib
import datetime
import json
import os
import re

import requests
from Bio import Entrez
from habanero import Crossref

Entrez.email = "talley.lambert+fpbase@gmail.com"
Entrez.api_key = os.getenv("NCBI_API_KEY", None)
email = Entrez.email
ID_CONVERT_URL = (
    "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=FPbase&email=" + email + "&ids=%s&format=json"
)


def pmc_converter(id, to="pmid"):
    """convert from Pubmed ID, Pubmed Central ID, or DOI, to anything else

    options for 'to':
        'pmid'  :   pubmed id
        'pmcid' :   pubmed central id
        'doi'   :   doi
        'all'   :   get full dict
    """
    many = False
    if isinstance(id, list | tuple | set):
        if len(id) > 200:
            out = []
            ids = [i for i in id if i]  # convert to list and remove null
            while len(ids):
                id_sub = ids[:200]
                out.extend(pmc_converter(id_sub, to=to))
                del ids[:200]
            return out
        else:
            id = ",".join(id)
            many = True
    record = json.loads(requests.get(ID_CONVERT_URL % id).content).get("records", None)
    try:
        if to == "all":
            if many:
                return record
            else:
                return record[0]
        else:
            if many:
                return [r.get(to) for r in record]
            else:
                return record[0].get(to, None)
    except (TypeError, IndexError):
        return None


def parse_crossref(doidict):
    if "message" in doidict:
        doidict = doidict["message"]
    out = {}
    ct = doidict.get("container-title", [])
    out = {
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


def crossref(doi):
    cr = Crossref(mailto="talley.lambert+fpbase@gmail.org")
    response = cr.works(ids=doi)
    # habanero returns a list if doi is a list of len > 1
    # otherwise a single dict
    if isinstance(doi, list | tuple | set) and len(doi) > 1:
        return {x.pop("doi"): x for x in (parse_crossref(i) for i in response)}
    else:
        return parse_crossref(response)


def doi2pmid(doi):
    pubmed_record = Entrez.read(Entrez.esearch(db="pubmed", term=doi))
    try:
        return pubmed_record.get("IdList")[0]
    except Exception:
        return None


def pmid2doi(pmid):
    pubmed_record = Entrez.read(Entrez.esummary(db="pubmed", id=pmid, retmode="xml"))
    try:
        return pubmed_record[0].get("DOI")
    except Exception:
        return None


def get_pmid_info(pmid):
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


def merge_info(dict1, dict2, exclude=()):
    """existings values in dict2 will overwrite dict1"""
    for key in dict1.keys():
        if key in dict2 and dict2[key] and key not in exclude:
            dict1[key] = dict2[key]
    return dict1


def doi_lookup(doi):
    info = crossref(doi)
    pmid = doi2pmid(doi)
    if pmid:
        try:
            pinfo = get_pmid_info(pmid)
            assert pinfo.pop("doi") == doi
            info = merge_info(pinfo, info)
            info["pmid"] = pmid
        except AssertionError:
            pass
    # get rid of empty values
    [info.pop(key) for key in list(info.keys()) if not info[key]]
    return info


def name_to_initials(name):
    return re.sub("([^A-Z-])", "", name)

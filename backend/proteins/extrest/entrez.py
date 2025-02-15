import json
import logging
import re
import time
import xml.etree.ElementTree as ET

from Bio import Entrez, SeqIO
from django.core.cache import cache

from references.helpers import pmid2doi

logger = logging.getLogger(__name__)

Entrez.email = "talley_lambert@hms.harvard.edu"


# genbank protein accession
gbprotrx = re.compile(r"^[A-Za-z]{3}\d{5}\.?\d?")
# genbank nucleotide accession
gbnucrx = re.compile(r"^([A-Za-z]{2}\d{6}\.?\d?)|([A-Za-z]{1}\d{5}\.?\d?)")
# refseq accession
refseqrx = re.compile("(NC|AC|NG|NT|NW|NZ|NM|NR|XM|XR|NP|AP|XP|YP|ZP)_[0-9]+")


def get_taxonomy_id(term, autochoose=False):
    record = Entrez.read(Entrez.esearch(db="taxonomy", term=term))
    if "ErrorList" in record and "PhraseNotFound" in record["ErrorList"]:
        spell_cor = get_entrez_spelling(term)
        if spell_cor:
            if autochoose:
                response = "y"
            else:
                response = input(f'By "{term}", did you mean "{spell_cor}"? (y/n): ')
            if response.lower() == "y":
                record = Entrez.read(Entrez.esearch(db="taxonomy", term=spell_cor))
    if record["Count"] == "1":
        return record["IdList"][0]
    elif int(record["Count"]) > 1:
        print("multiple records found:\n{}".format("\n".join(record["IdList"])))
    return None


def get_entrez_spelling(term, db="taxonomy"):
    record = Entrez.read(Entrez.espell(db=db, term=term))
    if record.get("CorrectedQuery"):
        return record.get("CorrectedQuery")
    else:
        return term


def get_ipgid_by_name(protein_name, give_options=True, recurse=True, autochoose=0):
    # autochoose value is the difference in protein counts required to autochoose
    # the most abundant protein count over the second-most abundant
    # example: autochoose=0 -> always choose the highest count, even in a tie
    #          autochoose=1 -> only chose if the highest is at least 1 greater...
    with Entrez.esearch(db="ipg", term=str(protein_name) + "[protein name]") as handle:
        record = Entrez.read(handle)
    if record["Count"] == "1":  # we got a unique hit
        return record["IdList"][0]
    elif record["Count"] == "0":
        # try without spaces
        if recurse:
            alternate_names = ["".join(protein_name.split()).lower()]
            alternate_names.append(alternate_names[0].replace("monomeric ", "m"))
            alternate_names.append(alternate_names[0].replace("-", ""))
            for name in set(alternate_names):
                if name != protein_name:
                    # print('Trying alternate name: {} -> {}'.format(protein_name, name))
                    uid = get_ipgid_by_name(name, recurse=False)
                    if uid:
                        return uid
            print(f"No results at IPG for: {protein_name}")
        return 0
    else:
        print("cowardly refusing to fetch squence with {} records at IPG".format(record["Count"]))
        if give_options:
            print("{:>4}{:<11}{:<30}{:<5}".format("", "ID", "NAME", "#PROT"))
            idlist = []
            for i, ipg_uid in enumerate(record["IdList"]):
                with Entrez.esummary(db="ipg", id=ipg_uid) as handle:
                    root = ET.fromstring(handle.read())
                docsum = root.find("DocumentSummarySet").find("DocumentSummary")
                prot_count = int(docsum.find("ProteinCount").text)
                name = docsum.find("Title").text
                print(f"{i + 1:>2}. {ipg_uid:<11}{name[:28]:<30}{prot_count:<5}")
                idlist.append((ipg_uid, prot_count))
            sorted_by_count = sorted(idlist, key=lambda tup: tup[1])
            sorted_by_count.reverse()
            diff = abs(sorted_by_count[0][1] - sorted_by_count[1][1])
            print("diff=", diff)
            if (autochoose >= 0) and (diff >= autochoose):
                outID = sorted_by_count[0][0]
                print(f"Autochoosing id {outID} with {sorted_by_count[0][1]} proteins:")
            else:
                rownum = input("Enter the row number you want to lookup, or press enter to cancel: ")
                try:
                    rownum = int(rownum) - 1
                    outID = idlist[rownum][0]
                except ValueError:
                    print("Exiting...")
            return outID
        return None


def fetch_ipg_sequence(protein_name=None, uid=None):
    """Retrieve protein sequence and IPG ID by name"""

    if not (protein_name or uid):
        raise ValueError("Must provide at least protein_name or uid")
    elif uid and not protein_name:
        ipg_uid = uid
    else:
        ipg_uid = get_ipgid_by_name(protein_name)

    try:
        with Entrez.esummary(db="ipg", id=ipg_uid) as handle:
            root = ET.fromstring(handle.read())
        docsum = root.find("DocumentSummarySet").find("DocumentSummary")
        prot_name = docsum.find("Title").text
    except AttributeError:
        return None

    print(f"Found protein with ID {ipg_uid}: {prot_name}")
    # prot_count = docsum.find('ProteinCount').text
    # assert prot_count == '1', 'Non-unique result returned'
    seq_len = docsum.find("Slen").text
    accession = docsum.find("Accession").text
    with Entrez.efetch(db="protein", id=accession, retmode="xml") as handle:
        record = Entrez.read(handle)
    assert len(record) == 1, "More than one record returned from protein database"
    record = record[0]
    prot_seq = record["GBSeq_sequence"].upper()
    assert len(prot_seq) == int(seq_len), (
        f"Protein database sequence different length {len(prot_seq)} than IPG database{int(seq_len)}"
    )
    return (ipg_uid, prot_seq)


def get_ipgid_from_gbid(gbid):
    with Entrez.esearch(db="ipg", term=gbid) as handle:
        record = Entrez.read(handle)
    if record["Count"] == "1":  # we got a unique hit
        return record["IdList"][0]


def check_accession_type(gbid):
    database = None
    if gbprotrx.match(gbid):  # it is a protein accession
        database = "protein"
    elif gbnucrx.match(gbid):
        database = "nuccore"
    elif gbid.upper().startswith(("WP_", "XP_", "YP_", "NP_", "AP_")):
        database = "protein"
    return database


def get_cached_gbseqs(gbids, max_age=60 * 60 * 24):
    gbseqs = cache.get("gbseqs", {})
    now = time.time()
    tofetch = [id for id in gbids if id not in gbseqs or gbseqs[id][1] - now >= max_age]
    gbseqs.update({k: (v, now) for k, v in fetch_gb_seqs(tofetch).items()})
    cache.set("gbseqs", gbseqs, 60 * 60 * 24)
    return gbseqs


def fetch_gb_seqs(gbids):
    """Retrieve protein sequence for multiple genbank IDs, (regardless of accession type)"""
    prots = []
    nucs = []
    for id in gbids:
        db = check_accession_type(id)
        if db == "protein":
            prots.append(id)
        elif db == "nuccore":
            nucs.append(id)
        else:
            prots.append(id)
            logger.error(f"Could not determine accession type for {id}")
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


def get_gb_seq(gbid):
    """Retrieve protein sequence for genbank ID, (regardless of accession type)"""
    database = check_accession_type(gbid)
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
        return parse_gbnuc_record(record).get("seq", None)
    return None


# ## MAIN FUNCTION ##
# use this to get lots of other info based on a genbank ID
def get_gb_info(gbid):
    """return dict of info for a given genbank id"""

    database = check_accession_type(gbid)
    if not database:
        return None

    with Entrez.efetch(db=database, id=gbid, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")

    if record:
        if database == "protein":
            D = parse_gbprot_record(record)
        elif database == "nuccore":
            D = parse_gbnuc_record(record)
        D["ipg_id"] = get_ipgid_from_gbid(gbid)
        if len(D["pmids"]):
            pmids = D.pop("pmids")
            D["refs"] = [{"pmid": pmid, "doi": pmid2doi(pmid)} for pmid in pmids]
        if D.get("organism"):
            taxid = get_taxonomy_id(D["organism"], True)
            if taxid:
                D["organism"] = int(taxid)
        if D.get("gb_prot"):
            from .uniprot import map_retrieve

            upids = map_retrieve(D.get("gb_prot"), source_fmt="EMBL").strip().split("\n")
            D["uniprots"] = [i for i in upids if i]
        return D
    else:
        return None


def parse_gbnuc_record(record):
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


def parse_gbprot_record(record):
    """
    returns (possible)
        desc: description
        seq: protein sequence
        pmids: [list of pmids]
        organism: name of organism
        gb_prot: protein accession
        gb_nuc:  nucleotide accession
        product: description of protein product
    """
    D = {}
    annotations = getattr(record, "annotations", None)
    features = getattr(record, "features", None)
    D["desc"] = getattr(record, "description", None)
    if hasattr(record, "_seq"):
        if hasattr(record._seq, "_data"):
            D["seq"] = record._seq._data

    if annotations:
        refs = annotations.get("references")
        if refs:
            D["pmids"] = [getattr(r, "pubmed_id", None) for r in refs if getattr(r, "pubmed_id", False)]
        D["organism"] = annotations.get("organism", None)
        acc = annotations.get("accessions", [])
        if len(acc):
            D["gb_prot"] = acc[0].split(".")[0]
        source = annotations.get("db_source")
        if source:
            D["gb_nuc"] = source.lstrip("accession ")  # noqa
    if features:
        prots = [f for f in features if f.type == "Protein"]
        if len(prots) == 1:
            prot = prots[0]
            quals = getattr(prot, "qualifiers", None)
            if quals:
                product = quals.get("product", [])
                if len(product):
                    D["product"] = product[0]
    return D


def get_otherid_from_ipgid(ipgid):
    with Entrez.esummary(db="ipg", id=ipgid, retmode="json") as handle:
        record = json.loads(handle.read())
        try:
            return record["result"][ipgid]["accession"]
        except Exception:
            return None

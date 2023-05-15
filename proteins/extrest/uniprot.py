from xml.dom import minidom

import requests

from proteins.validators import validate_uniprot

BASE = "http://www.uniprot.org"
KB_ENDPOINT = "/uniprot/"
TOOL_ENDPOINT = "/uploadlists/"


def get_uniprot_info(unip_id):
    """Get sequence, accessions, genbank matches, and references from Uniprot"""
    try:
        validate_uniprot(unip_id)
    except Exception:
        print(f"ERROR, invalid uniprot ID: {unip_id}")
        return {}
    D = parse_uniprot_record(get_uniprot_record(unip_id))
    D["uniprot"] = unip_id
    gbids = map_retrieve(unip_id, target_fmt="EMBL").strip().split("\n")
    D["genbank"] = []
    for id in gbids:
        from .entrez import get_gb_seq

        if get_gb_seq(id) == D["seq"]:
            D["genbank"].append(id.split(".")[0])
    return D


def get_uniprot_record(id):
    response = requests.get(BASE + KB_ENDPOINT + str(id) + ".xml")
    if response.status_code == 200:
        return response.content
    else:
        return None


def parse_uniprot_record(content):
    xmldoc = minidom.parseString(content)
    D = {"refs": [], "seq": None, "all_uniprots": [], "pdbs": []}

    getE = xmldoc.getElementsByTagName
    D["all_uniprots"] = [item.firstChild.data for item in getE("accession")]
    for dbR in getE("dbReference"):
        if dbR.getAttribute("type") == "DOI":
            D["refs"].append(dbR.getAttribute("id"))
        elif dbR.getAttribute("type") == "PDB":
            D["pdbs"].append(dbR.getAttribute("id"))
    D["seq"] = getE("sequence")[0].firstChild.data.replace("\n", "")

    return D


def map_retrieve(ids2map, source_fmt="ACC+ID", target_fmt="ACC", output_fmt="list"):
    """Map database identifiers from/to UniProt accessions.
    The mapping is achieved using the RESTful mapping service provided by
    UniProt. While a great many identifiers can be mapped the documentation
    has to be consulted to check which options there are and what the database
    codes are. Mapping UniProt to UniProt effectlvely allows batch retrieval
    of entries.
    Args:
    ids2map (list or string): identifiers to be mapped
    source_fmt (str, optional): format of identifiers to be mapped.
    Defaults to ACC+ID, which are UniProt accessions or IDs.
    target_fmt (str, optional): desired identifier format. Defaults
    to ACC, which is UniProt accessions.
    output_fmt (str, optional): return format of data. Defaults to list.
    Returns:
    mapped identifiers (str)
    """
    if hasattr(ids2map, "pop"):
        ids2map = " ".join(ids2map)
    payload = {
        "from": source_fmt,
        "to": target_fmt,
        "format": output_fmt,
        "query": ids2map,
    }
    response = requests.get(BASE + TOOL_ENDPOINT, params=payload)
    if response.ok:
        return response.text
    else:
        response.raise_for_status()

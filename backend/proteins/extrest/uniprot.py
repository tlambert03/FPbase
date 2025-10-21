from xml.dom import minidom

from external_apis import sequences
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
    try:
        xml_text = sequences.fetch_uniprot_xml(id)
        return xml_text.encode("utf-8")  # Return bytes to match original behavior
    except Exception:
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
        ids2map = list(ids2map)
    else:
        ids2map = [ids2map]
    return sequences.map_uniprot_ids(ids2map, from_db=source_fmt, to_db=target_fmt)

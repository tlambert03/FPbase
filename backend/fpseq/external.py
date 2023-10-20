import json

import requests


def genbank_seq(accession):
    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        f"efetch.fcgi?db=protein&id={accession}&rettype=fasta&"
        "retmode=text"
    )
    if response:
        return "".join([x.decode() for x in response.content.splitlines()[1:]])
    return None


def uniprot_seq(accession):
    response = requests.get(f"https://www.uniprot.org/uniprot/{accession}.fasta")
    if response:
        return "".join([x.decode() for x in response.content.splitlines()[1:]])
    return None


def pdb_seq(accession):
    response = requests.get(f"http://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{accession}")
    if response:
        j = json.loads(response.text)
        try:
            j = j[accession.lower()]
        except ValueError:
            j = j[accession]
        return j[0]["sequence"]
    return None

import requests
import json


def genbank_seq(accession):
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
                            'efetch.fcgi?db=protein&id={}&rettype=fasta&'
                            'retmode=text'.format(accession))
    if response:
        return "".join([x.decode() for x in response.content.splitlines()[1:]])
    return None


def uniprot_seq(accession):
    response = requests.get('https://www.uniprot.org/uniprot/{}.fasta'.format(accession))
    if response:
        return "".join([x.decode() for x in response.content.splitlines()[1:]])
    return None


def pdb_seq(accession):
    response = requests.get('http://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{}'.format(accession))
    if response:
        j = json.loads(response.text)
        try:
            j = j[accession.lower()]
        except ValueError:
            j = j[accession]
        return j[0]['sequence']
    return None

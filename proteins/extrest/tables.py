from bs4 import BeautifulSoup
import requests
import tablib


def parensplit(text):
    v = text.strip(')').split('(')
    if len(v) == 2:
        val, unit = v
    else:
        val = " ".join(v)
        unit = ''
    return (val.strip(), unit)


def table2dataset(table):
    data = tablib.Dataset()
    data.headers = [head.text for head in table.find('thead').find_all('th')]
    for row in table.find('tbody').find_all('tr'):
        data.append(tuple([td.text.strip() for td in row.find_all('td')]))
    return data


def fetch_pmc_content(pmcid):
    URI = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id='
    return requests.get(URI + pmcid)


def fetch_doi_content(doi):
    return requests.get('https://doi.org/' + doi)


def response2table2(response):
    soup = BeautifulSoup(response.content, 'lxml')
    tables = soup.find_all('table')
    return [table2dataset(table) for table in tables]


def pmcid2tables(pmcid):
    response = fetch_pmc_content(pmcid)
    if not response.status_code == 200:
        print('Bad response: {}'.format(response.status_code))
        return None
    return response2table2(response)


def doi2tables(doi):
    response = fetch_doi_content(doi)
    if not response.status_code == 200:
        print('Bad response: {}'.format(response.status_code))
        return None
    return response2table2(response)


dois = set([p.primary_reference.doi for p in Protein.objects.all() if p.primary_reference])
pmcids = [p for p in pmc_converter(dois, to='pmcid') if p]


class HMSProxy(EzProxy):
    def __init__(self, barcode, last_name):
        super().__init__("ezp-prod1.hul.harvard.edu")
        self.login(barcode, last_name)



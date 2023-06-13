import re

import pandas as pd
import requests
import tablib  # TODO: convert to pandas
from bs4 import BeautifulSoup


def interpret_heading(head_str):
    HEADING_RX = {
        "protein": r"protein",
        "abs_max": r"absorb",
        "ex_max": r"(excitation|λex)",
        "em_max": r"(emission|λem)",
        "stability": r"(bleach|photostab)",
        "maturation": r"matur",
        "brightness": r"bright",
        "lifetime": r"lifetime",
        "ext_coeff": r"(M[--]1\s*cm[--]1|ɛ|ε|extinction|^ec\s)",
        "QY": r"(^qy|quantum|ϕ)",
        "pka": r"pka",
        "agg": r"oligomer",
    }

    matches = []
    for key, rgx in HEADING_RX.items():
        if re.search(rgx, head_str, re.IGNORECASE):
            matches.append(key)
    if not matches:
        return None
    elif len(matches) > 1:
        print("MULTIPLE MATCHES!!!")
    return matches[0]


def parensplit(text):
    v = text.strip(")").split("(")
    if len(v) == 2:
        val, unit = v
    else:
        val = " ".join(v)
        unit = ""
    return (val.strip(), unit)


def first_row(table):
    return [th.get_text() for th in table.find("tr").find_all("th")]


def table2dataset(table):
    if isinstance(table, str):
        table = BeautifulSoup(table, "lxml").find_all("table")[0]
    data = tablib.Dataset()
    # data.headers = [head.text for head in table.find('thead').find_all('th')]
    headings = first_row(table)
    data.headers = headings
    for row in table.find("tbody").find_all("tr"):
        data.append(tuple([td.text.strip() for td in row.find_all("td")]))
    return data


def fetch_pmc_content(pmcid):
    URI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id="
    return requests.get(URI + pmcid)


def fetch_doi_content(doi):
    return requests.get("https://doi.org/" + doi)


def text2tables(text):
    soup = BeautifulSoup(text, "lxml")
    tables = soup.find_all("table")
    print(f"found {len(tables)} tables")
    return tables


def response2table2(response):
    tables = text2tables(response.content)
    return [table2dataset(table) for table in tables]


def pmcid2tables(pmcid):
    response = fetch_pmc_content(pmcid)
    if not response.status_code == 200:
        print(f"Bad response: {response.status_code}")
        return None
    return response2table2(response)


def doi2tables(doi):
    response = fetch_doi_content(doi)
    if not response.status_code == 200:
        print(f"Bad response: {response.status_code}")
        return None
    return response2table2(response)


class HTMLTableParser:
    def parse_text(self, text):
        soup = BeautifulSoup(text, "lxml")
        return [self.parse_html_table(table) for table in soup.find_all("table")]

    def parse_url(self, url):
        response = requests.get(url)
        return self.parse_text(response.text)

    def parse_html_table(self, table):
        n_columns = 0
        n_rows = 0
        column_names = []

        # Find number of rows and columns
        # we also find the column titles if we can
        for row in table.find_all("tr"):
            # Determine the number of rows in the table
            td_tags = row.find_all("td")
            if len(td_tags) > 0:
                n_rows += 1
                if n_columns == 0:
                    # Set the number of columns for our table
                    n_columns = len(td_tags)

            # Handle column names if we find them
            th_tags = row.find_all("th")
            if len(th_tags) > 0 and len(column_names) == 0:
                for th in th_tags:
                    column_names.append(th.get_text())

        # Safeguard on Column Titles
        if len(column_names) > 0 and len(column_names) != n_columns:
            raise Exception("Column titles do not match the number of columns")

        columns = column_names if len(column_names) > 0 else range(0, n_columns)
        df = pd.DataFrame(columns=columns, index=range(0, n_rows))
        row_marker = 0
        for row in table.find_all("tr"):
            column_marker = 0
            columns = row.find_all("td")
            for column in columns:
                df.iat[row_marker, column_marker] = column.get_text()
                column_marker += 1
            if len(columns) > 0:
                row_marker += 1

        # Convert to float if possible
        for col in df:
            try:
                df[col] = df[col].astype(float)
            except ValueError:
                pass

        return df


# dois = set([p.primary_reference.doi for p in Protein.objects.all() if p.primary_reference])
# pmcids = [p for p in pmc_converter(dois, to='pmcid') if p]


# class HMSProxy(EzProxy):
#     def __init__(self, barcode, last_name):
#         super().__init__("ezp-prod1.hul.harvard.edu")
#         self.login(barcode, last_name)

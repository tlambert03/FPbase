import requests
import json
import datetime
from fpbase import __version__
import re

from Bio import Entrez
from crossref.restful import Works, Etiquette
Entrez.email = 'talley.lambert+fpbase@gmail.com'
email = Entrez.email
my_etiquette = Etiquette('FPbase', __version__, 'http://fpbase.org', email)
ID_CONVERT_URL = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=FPbase&email='+email+'&ids=%s&format=json'


def pmc_converter(id, to='pmid',):
    """convert from Pubmed ID, Pubmed Central ID, or DOI, to anything else

    options for 'to':
        'pmid'  :   pubmed id
        'pmcid' :   pubmed central id
        'doi'   :   doi
        'all'   :   get full dict
    """
    many = False
    if isinstance(id, (list, tuple, set)):
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
    record = json.loads(requests.get(ID_CONVERT_URL % id).content).get('records', None)
    try:
        if to == 'all':
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


def get_doi_info(doi):
    works = Works(etiquette=my_etiquette)
    query = works.doi(doi)
    try:
        journal = query['container-title'][0]
    except Exception:
        journal = None
    try:
        title = query.get('title', None)
        if title:
            title = title[0]
    except Exception:
        title = None
    try:
        journal_iss = query.get('journal-issue', {})
        dateparts = journal_iss.get('published-online', {})
        if not dateparts:
            dateparts = journal_iss.get('published-print', {})

        try:
            dateparts = dateparts.get('date-parts', None)[0]
        except Exception:
            pass

        if dateparts:
            pubyear = dateparts[0]
            if len(dateparts) == 1:
                dateparts.append(1)
            if len(dateparts) == 2:
                dateparts.append(1)
            dateparts = datetime.date(*dateparts)
        else:
            issued = query.get('issued', [])
            if issued:
                pubyear = issued['date-parts'][0][0]
                dp = issued['date-parts'][0]
                if len(dp) == 1:
                    dp.append(1)
                if len(dp) == 2:
                    dp.append(1)
                dateparts = datetime.date(*dp)
    except Exception as e:
        pubyear = None
        dateparts = None
        print(e)
    return {
        'title': title,
        'journal': journal,
        'pages': query.get('page', None),
        'volume': query.get('volume', None),
        'issue': query.get('issue', None),
        'authors': query.get('author', None),
        'year': pubyear,
        'date': dateparts
    }


def doi2pmid(doi):
    pubmed_record = Entrez.read(Entrez.esearch(db='pubmed', term=doi))
    try:
        return pubmed_record.get('IdList')[0]
    except Exception:
        return None


def pmid2doi(pmid):
    pubmed_record = Entrez.read(Entrez.esummary(db='pubmed', id=pmid, retmode='xml'))
    try:
        return pubmed_record[0].get('DOI')
    except Exception:
        return None


def get_pmid_info(pmid):
    pubmed_record = Entrez.read(Entrez.esummary(db='pubmed', id=pmid, retmode='xml'))
    if len(pubmed_record):
        pubmed_record = pubmed_record[0]
        date = None
        try:
            date = datetime.datetime.strptime(pubmed_record['PubDate'], "%Y %b %d").date()
        except Exception:
            try:
                date = datetime.datetime.strptime(pubmed_record['EPubDate'], "%Y %b %d").date()
            except Exception:
                pass
        return {
            'doi': pubmed_record.get('DOI', None),
            'title': pubmed_record.get('Title', ''),
            'journal': pubmed_record.get('Source', ''),
            'pages': pubmed_record.get('Pages', ''),
            'volume': pubmed_record.get('Volume', ''),
            'issue': pubmed_record.get('Issue', ''),
            'year': pubmed_record['PubDate'].split()[0],
            'authors': pubmed_record['AuthorList'],
            'date': date,
        }


def merge_info(dict1, dict2, exclude=[]):
    """ existings values in dict2 will overwrite dict1 """
    for key in dict1.keys():
        if key in dict2 and dict2[key]:
            if key not in exclude:
                dict1[key] = dict2[key]
    return dict1


def doi_lookup(doi):
    info = get_doi_info(doi)
    pmid = doi2pmid(doi)
    if pmid:
        try:
            pinfo = get_pmid_info(pmid)
            assert pinfo.pop('doi') == doi
            info = merge_info(pinfo, info)
            info['pmid'] = pmid
        except AssertionError:
            pass
    # get rid of empty values
    [info.pop(key) for key in list(info.keys()) if not info[key]]
    return info


def name_to_initials(name):
    return re.sub('([^A-Z-])', '', name)

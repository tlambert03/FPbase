import requests
import json
from crossref.restful import Works, Etiquette
from fpbase import __version__
import re

from Bio import Entrez
Entrez.email = 'talley.lambert+fpbase@gmail.com'
email = Entrez.email
my_etiquette = Etiquette('FPbase', __version__, 'http://fpbase.org', email)
ID_CONVERT_URL = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=FPbase&email='+email+'&ids=%s&format=json'


def get_pmid_from_doi(doi):
	record = json.loads(requests.get(ID_CONVERT_URL % doi).content).get('records', None)
	try:
		return record[0].get('pmid', None)
	except IndexError:
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
		pubyear = query.get('issued', None)
		if pubyear:
			pubyear = pubyear['date-parts'][0][0]
	except Exception as e:
		pubyear = None
		print(e)
	return {
		'title': title,
		'journal': journal,
		'pages': query.get('page', None),
		'volume': query.get('volume', None),
		'issue': query.get('issue', None),
		'authors': query.get('author', None),
		'year': pubyear,
	}


def get_doi_from_pmid(pmid):
	pubmed_record = Entrez.read(Entrez.esummary(db='pubmed', id=pmid, retmode='xml'))
	try:
		return pubmed_record[0].get('DOI')
	except Exception:
		return None


def get_pmid_info(pmid):
	pubmed_record = Entrez.read(Entrez.esummary(db='pubmed', id=pmid, retmode='xml'))
	if len(pubmed_record):
		pubmed_record = pubmed_record[0]
		return {
			'doi': pubmed_record.get('DOI', None),
			'title': pubmed_record.get('Title', ''),
			'journal': pubmed_record.get('Source', ''),
			'pages': pubmed_record.get('Pages', ''),
			'volume': pubmed_record.get('Volume', ''),
			'issue': pubmed_record.get('Issue', ''),
			'year': pubmed_record['PubDate'].split()[0],
			'authors': pubmed_record['AuthorList']
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
	pmid = get_pmid_from_doi(doi)
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

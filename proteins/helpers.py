import xml.etree.ElementTree as ET
from metapub import CrossRef, PubMedFetcher
from Bio import Entrez
import re
CR = CrossRef()
PMF = PubMedFetcher()
# Entrez.parse doesn't seem to work with ipg results

Entrez.email = "talley_lambert@hms.harvard.edu"


def wave_to_hex(wavelength, gamma=1):
    '''This converts a given wavelength into an approximate RGB value.
    The given wavelength is in nanometers.
    The range of wavelength is 380 nm through 750 nm.

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''

    wavelength = float(wavelength)
    if 520 <= wavelength:
        #pass
        wavelength += 40

    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 850:
        attenuation = 0.3 + 0.7 * (770 - wavelength) / (770 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    R *= 255
    G *= 255
    B *= 255
    return '#%02x%02x%02x' % (int(R), int(G), int(B))


def get_color_group(ex_max, em_max):
    if (em_max - ex_max) > 90:
        return "Long Stokes Shift", "#80A0FF"
    if ex_max < 380:
        return "UV", "#C080FF"
    if ex_max < 421:
        return "Blue", "#8080FF"
    if ex_max < 473:
        return "Cyan", "#80FFFF"
    if ex_max < 505:
        return "Green", "#80FF80"
    if ex_max < 515:
        return "Green/Yellow", "#CCFF80"
    if ex_max < 531:
        return "Yellow", "#FFFF80"
    if ex_max < 555:
        return "Orange", "#FFC080"
    if ex_max < 600:
        return "Red", "#FFA080"
    if ex_max < 631:
        return "Far Red", "#FF8080"
    if ex_max < 800:
        return "Near IR", "#B09090"


def get_taxonomy_id(term, autochoose=False):
    record = Entrez.read(Entrez.esearch(db='taxonomy', term=term))
    if 'ErrorList' in record and 'PhraseNotFound' in record['ErrorList']:
        spell_cor = get_entrez_spelling(term)
        if spell_cor:
            if autochoose:
                response = 'y'
            else:
                response = input('By "{}", did you mean "{}"? (y/n): '.format(term, spell_cor))
            if response.lower() == 'y':
                record = Entrez.read(Entrez.esearch(db='taxonomy', term=spell_cor))
    if record['Count'] == '1':
        return record['IdList'][0]
    elif int(record['Count']) > 1:
        print('multiple records found:\n{}'.format('\n'.join(record['IdList'])))
    return None


def get_entrez_spelling(term, db='taxonomy'):
    record = Entrez.read(Entrez.espell(db=db, term=term))
    if record.get('CorrectedQuery'):
        return record.get('CorrectedQuery')
    else:
        return term


def get_ipgid_by_name(protein_name, give_options=True, recurse=True, autochoose=0):
    # autochoose value is the difference in protein counts required to autochoose
    # the most abundant protein count over the second-most abundant
    # example: autochoose=0 -> always choose the highest count, even in a tie
    #          autochoose=1 -> only chose if the highest is at least 1 greater...
    with Entrez.esearch(db='ipg', term=str(protein_name) + '[protein name]') as handle:
        record = Entrez.read(handle)
    if record['Count'] == '1':  # we got a unique hit
        return record['IdList'][0]
    elif record['Count'] == '0':
        # try without spaces
        if recurse:
            alternate_names = ["".join(protein_name.split()).lower()]
            alternate_names.append(alternate_names[0].replace('monomeric ', 'm'))
            alternate_names.append(alternate_names[0].replace('-', ''))
            for name in set(alternate_names):
                if name != protein_name:
                    # print('Trying alternate name: {} -> {}'.format(protein_name, name))
                    uid = get_ipgid_by_name(name, recurse=False)
                    if uid:
                        return uid
            print('No results at IPG for: {}'.format(protein_name))
        return 0
    else:
        print('cowardly refusing to fetch squence with {} records at IPG'.format(record['Count']))
        if give_options:
            print('{:>4}{:<11}{:<30}{:<5}'.format("", "ID", "NAME", "#PROT"))
            idlist = []
            for i, ipg_uid in enumerate(record['IdList']):
                with Entrez.esummary(db='ipg', id=ipg_uid) as handle:
                    root = ET.fromstring(handle.read())
                docsum = root.find('DocumentSummarySet').find('DocumentSummary')
                prot_count = int(docsum.find('ProteinCount').text)
                name = docsum.find('Title').text
                print('{:>2}. {:<11}{:<30}{:<5}'.format(i + 1, ipg_uid, name[:28], prot_count))
                idlist.append((ipg_uid, prot_count))
            sorted_by_count = sorted(idlist, key=lambda tup: tup[1])
            sorted_by_count.reverse()
            diff = abs(sorted_by_count[0][1] - sorted_by_count[1][1])
            print('diff=', diff)
            if (autochoose >= 0) and (diff >= autochoose):
                outID = sorted_by_count[0][0]
                print('Autochoosing id {} with {} proteins:'.format(outID, sorted_by_count[0][1]))
            else:
                rownum = input('Enter the row number you want to lookup, or press enter to cancel: ')
                try:
                    rownum = int(rownum) - 1
                    outID = idlist[rownum][0]
                except ValueError:
                    print('Exiting...')
            return outID
        return None


def fetch_ipg_sequence(protein_name=None, uid=None):
    """Retrieve protein sequence and IPG ID by name"""

    if not (protein_name or uid):
        raise ValueError('Must provide at least protein_name or uid')
    elif uid and not protein_name:
        ipg_uid = uid
    else:
        ipg_uid = get_ipgid_by_name(protein_name)

    try:
        with Entrez.esummary(db='ipg', id=ipg_uid) as handle:
            root = ET.fromstring(handle.read())
        docsum = root.find('DocumentSummarySet').find('DocumentSummary')
        prot_name = docsum.find('Title').text
    except AttributeError:
        return None

    print("Found protein with ID {}: {}".format(ipg_uid, prot_name))
    # prot_count = docsum.find('ProteinCount').text
    # assert prot_count == '1', 'Non-unique result returned'
    seq_len = docsum.find('Slen').text
    accession = docsum.find('Accession').text
    with Entrez.efetch(db='protein', id=accession, retmode='xml') as handle:
        record = Entrez.read(handle)
    assert len(record) == 1, 'More than one record returned from protein database'
    record = record[0]
    prot_seq = record['GBSeq_sequence'].upper()
    assert len(prot_seq) == int(seq_len), 'Protein database sequence different length {} than IPG database{}'.format(len(prot_seq), int(seq_len))
    return (ipg_uid, prot_seq)


def mless(name):
    if re.search('^m[A-Z]', name):
        return name.lstrip('m')
    if name.startswith('monomeric'):
        name = name.lstrip('monomeric')
    if name.startswith('Monomeric'):
        name = name.lstrip('Monomeric')
    return name.lstrip(' ')


def get_base_name(name):
    '''return core name of protein, stripping prefixes like "m" or "Tag"'''

    # remove PA/(Pa), PS, PC, from beginning
    if name.startswith(('PA', 'Pa', 'PS', 'Ps', 'PC', 'pc', 'rs')):
        name = name[2:]

    if re.match('LSS', name):
        name = name[3:].lstrip('-')

    # remove m (if next letter is caps) or monomeric
    if re.match('m[A-Z]', name):
        name = name[1:]

    # get rid of Td or td
    if re.match('[Tt][Dd][A-Z]', name):
        name = name[2:]

    name = name.lstrip('Monomeric')
    name = name.lstrip('Tag')

    # remove E at beginning (if second letter is caps)
    if re.match('E[A-Z]', name):
        name = name[1:]
    # remove S at beginning (if second letter is caps)
    if re.match('S[A-Z]', name):
        name = name[1:]
    # remove T- at beginning (if second letter is caps)
    if re.match('T-', name):
        name = name[2:]

    name = name.lstrip('-').lstrip(' ')

    return name


import xml.etree.ElementTree as ET
import pandas as pd
import traceback

from proteins.models import Protein, State, StateTransition, BleachMeasurement
from references.models import Reference
from fpbase.users.models import User
from django.template.defaultfilters import slugify
from Bio import Entrez

from metapub import CrossRef, PubMedFetcher
CR = CrossRef()
PMF = PubMedFetcher()
# Entrez.parse doesn't seem to work with ipg results


Entrez.email = "talley_lambert@hms.harvard.edu"

############################################
#       Importing Tools
############################################

SUPERUSER = User.objects.filter(is_superuser=True).first()


def require_superuser(func):
    def wrapper(*args, **kwargs):
        if not SUPERUSER:
            raise ValueError('No SUPERUSER in database... cannot run function: {}'.format(func.__name__))
        return func(*args, **kwargs)
    return wrapper


aggLookup = {
    'Tetramer': 't',
    'Monomer': 'm',
    'Dimer': 'd'
}


def get_nonan(obj, item):
    val = obj.get(item, None)
    if val is not None:
        if pd.isna(val):
            return None
    return val


@require_superuser
def add_ref_to_prot(protein, doi, showexisting=False):
    ref, created = Reference.objects.get_or_create(doi=doi)
    rf = 1 if created else 0
    protein.references.add(ref)
    if not protein.primary_reference:
        protein.primary_reference = ref
    protein.save()
    return rf


@require_superuser
def importCSV(file=None):
    '''
    mainly intended as an import function for the proteins in the NIC table
    '''
    if file is None:
        url = 'https://raw.githubusercontent.com/FPvisualization/fpvisualization.github.io/master/FPs.csv'
        df = pd.read_csv(url)
    else:
        df = pd.read_csv(file)

    ps = 0
    st = 0
    rf = 0
    for i, prot in df.iterrows():
        if not Protein.objects.filter(slug=slugify(prot.Name)).exists():
            print("importing {}...".format(prot.Name))
            p = Protein(
                name=prot.Name,
                created_by=SUPERUSER,
                agg=get_nonan(prot, 'agg'),
                switch_type=Protein.BASIC,
            )
            p.save()
            ps += 1

            s = State(
                # name        = 'default',  # now the default state name
                ex_max=get_nonan(prot, 'lambda_ex'),
                em_max=get_nonan(prot, 'lambda_em'),
                ext_coeff=get_nonan(prot, 'E'),
                qy=get_nonan(prot, 'QY'),
                pka=get_nonan(prot, 'pka'),
                maturation=get_nonan(prot, 'mature'),
                lifetime=get_nonan(prot, 'lifetime'),
                protein=p,
                created_by=SUPERUSER,
            )
            s.save()
            st += 1

            # add it as default state
            p.default_state = s
            p.save()

            # add bleach numbers
            if get_nonan(prot, 'bleach'):
                BleachMeasurement.objects.update_or_create(
                    rate=get_nonan(prot, 'bleach'),
                    state=s,
                )
        else:
            print("{} already in database...".format(prot.Name))
            p = Protein.objects.get(slug=slugify(prot.Name))

        try:
            if pd.isna(prot.DOI) or not prot.DOI.startswith('10'):
                continue
            for doi in prot.DOI.split(" "):
                if not (doi and doi.startswith('10')):
                    continue
                rf += add_ref_to_prot(p, doi)
        except Exception as e:
            # traceback.print_exc()
            print("error importing reference: {}".format(e))

    print("{} Proteins, {} States, and {} References imported".format(ps, st, rf))


def linkstates(df):
    linksurl = 'https://raw.githubusercontent.com/FPvisualization/fpvisualization.github.io/master/links.csv'
    linksdf = pd.read_csv(linksurl)
    q = df.set_index('UID').to_dict()

    for i, link in linksdf.iterrows():
        p = Protein.objects.get(slug=slugify(q['Name'][link.state1]))
        fromState = p.states.get(name=q['state'][link.state1])
        toState = p.states.get(name=q['state'][link.state2])

        t, created = StateTransition.objects.get_or_create(
            protein=p,
            trans_wave=int(link.lambda_sw),
            from_state=fromState,
            to_state=toState,
        )
        if created:
            print('created: {}'.format(t))


@require_superuser
def importPSFPs(file=None):
    if file is None:
        url = 'https://raw.githubusercontent.com/kthorn/FPvisualization/master/PSFPs.csv'
        df = pd.read_csv(url)
    else:
        df = pd.read_csv(file)

    ps = 0
    st = 0
    rf = 0
    for i, prot in df.iterrows():
        try:

            p, created = Protein.objects.get_or_create(
                slug=slugify(prot.Name),
                defaults={
                    'name': prot.Name,
                    'created_by': SUPERUSER,
                    'switch_type': prot.type,
                    'agg': aggLookup[prot.Aggregation] if get_nonan(prot, 'Aggregation') else None
                }
            )
            if created:
                print("PROTEIN CREATED: {}".format(prot.Name))
                ps += 1
            else:
                print("Protein Found  : {}".format(prot.Name))

            if not pd.isna(prot.DOI) and prot.DOI.startswith('10'):
                try:
                    rf += add_ref_to_prot(p, prot.DOI)
                except Exception as e:
                    # traceback.print_exc()
                    print('Error importing reference: {}'.format(e))

            state, created = p.states.get_or_create(name=prot.state,
                                   defaults={
                                       'name': prot.state,
                                       'is_dark': 'True' if 'off' in prot.state.lower() else False,
                                       'ex_max': get_nonan(prot, 'lambda_ex'),
                                       'em_max': get_nonan(prot, 'lambda_em'),
                                       'ext_coeff': get_nonan(prot, 'E'),
                                       'qy': get_nonan(prot, 'QY'),
                                       'pka': get_nonan(prot, 'pka'),
                                       'maturation': get_nonan(prot, 'mature'),
                                       'lifetime': get_nonan(prot, 'lifetime'),
                                       'protein': p,
                                       'created_by': SUPERUSER,
                                   })

            # needs work... would like ON state to show up in tables...
            # but only if it is a DARK -> ON protein
            if p.switch_type == 'ps' or p.switch_type == 'pa':
                if not state.is_dark:
                    p.default_state = state
                    p.save()
            elif prot.get('initialState', False):
                p.default_state = state
                p.save()

            if get_nonan(prot, 'bleach'):
                BleachMeasurement.objects.update_or_create(
                    rate=get_nonan(prot, 'bleach'),
                    state=state,
                )
            if created:
                print("STATE CREATED  : {}".format(prot.state))
                st += 1
            else:
                print("State Found    : {}".format(prot.state))

        except Exception as e:
            traceback.print_exc()
            print('failed to import {}: {}'.format(prot.Name, e))

    print("{} Proteins, {} States, and {} References imported".format(ps, st, rf))
    linkstates(df)


def importSeqs(file=None):
    if file is None:
        import os
        basedir = os.path.dirname(os.path.dirname(__file__))
        url = os.path.join(basedir, '_data/FPseqs.csv')
        df = pd.read_csv(url)
    else:
        df = pd.read_csv(file)

    rf = 0
    sq = 0
    for i, prot in df.iterrows():
        if Protein.objects.filter(name__icontains=prot.Name).count() == 1:
            p = Protein.objects.get(name__icontains=prot.Name)
            if p.seq is None:
                p.seq = prot.AminoAcidSequence
                p.save()
                print("Added sequence to {}".format(prot.Name))
                sq += 1
            else:
                seq = prot.AminoAcidSequence.upper()
                seq = "".join(seq.split())
                if p.seq.upper() != seq:
                    print("Non-matching sequence found for {}!".format(prot.Name))
            try:
                if 'dx.doi' in prot.Source:
                    doi = prot.Source.strip('http://dx.doi.org/')
                    rf += add_ref_to_prot(p, doi)
            except Exception as e:
                # traceback.print_exc()
                print('Error importing reference: {}'.format(e))

    print("{} Sequences added; {} References imported".format(sq, rf))


def importSpectra(file=None):
    if file is None:
        import os
        basedir = os.path.dirname(os.path.dirname(__file__))
        print(basedir)
        url = os.path.join(basedir, '_data/FLUOR.csv')
        df = pd.read_csv(url)
    else:
        df = pd.read_csv(file)

    sp = 0
    for i, prot in df.iterrows():
        if Protein.objects.filter(name=prot.fluor_name).count() == 1:
            p = Protein.objects.get(name=prot.fluor_name)
            if not p.states.count() == 1:
                print('avoiding protein with multiple states: {}'.format(p))
                continue
            try:
                D = p.default_state
                if not D.ex_spectra:
                    D.ex_spectra = prot.ex_spectra
                    sp += 1
                if not D.em_spectra:
                    D.em_spectra = prot.em_spectra
                    sp += 1
                D.save()
            except Exception as e:
                print("failed to import spectrum for {}".format(prot.fluor_name))
                print(e)

    print("Imported {} spectra".format(sp))


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


@require_superuser
def reload_all():
    importCSV()
    importPSFPs()
    importSeqs()
    importSpectra()
    for P in Protein.objects.all():
        Q = fetch_ipg_sequence(P.name)
        if Q:
            P.seq = Q[1]
            P.ipg_id = Q[0]
            P.save()


#########################
#  Conversion
#########################

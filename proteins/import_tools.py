import pandas as pd
import traceback
import os
import tablib
import json
from proteins.models import Protein, State, StateTransition, BleachMeasurement, Organism
from references.models import Reference
from fpbase.users.models import User
from django.template.defaultfilters import slugify
from Bio import Entrez
from .helpers import fetch_ipg_sequence
from .validators import protein_sequence_validator
from proteins import forms

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
    'Dimer': 'd',
    'Dimer/Trimer': 'd',
}

FIELD_LOOKUP = {
    'chromophore_class': (Protein, 'chromophore'),
    'protein_name': (Protein, 'name'),
    'species': (Organism, 'scientific_name'),
    'excitation_maxima': (State, 'ex_max'),
    'emission_maxima': (State, 'em_max'),
    'quantum_yield': (State, 'qy'),
    'extinction_coefficient': (State, 'ext_coeff'),
    'oligomerization_state': ('protein', 'agg'),
    'pka': ('protein', 'pka'),
    'maturation': ('protein', 'maturation'),
    'bleaching': ('protein', 'rate'),
    'fluorescence_lifetime': ('protein', 'lifetime'),
    'doi': ('protein', 'doi'),
    'amino_acid_sequence': ('protein', 'seq'),
}


def get_nonan(obj, item):
    val = obj.get(item, None)
    if val is not None:
        if pd.isna(val):
            return None
    return val


@require_superuser
def add_ref_to_prot(protein, doi, showexisting=False):
    ref, created = Reference.objects.get_or_create(doi=doi.strip())
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
        basedir = os.path.dirname(os.path.dirname(__file__))
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


def import_organisms():
    with open('_data/species.json', 'r') as f:
        D = json.load(f)
    for k in D.keys():
        o = Organism(id=k)
        o.save()


ORGLOOKUP = {'Acanthastrea sp. ': None,
 'Acropora aculeus': 287157,
 'Acropora digitifera': 70779,
 'Acropora eurostoma': 526283,
 'Acropora hyacinthus': 55974,
 'Acropora millepora': 45264,
 'Acropora nobilis': 70781,
 'Acropora pulchra': 140239,
 'Acropora sp.': None,
 'Acropora tenuis': 70783,
 'Actinia equina': 6106,
 'Actinia equina ': 6106,
 'Aequorea Victoria': 6100,
 'Aequorea Victoria ': 6100,
 'Aequorea Victoria  ': 6100,
 'Aequorea coerulescens': 210840,
 'Aequorea coerulescens ': 210840,
 'Aequorea macrodactyla': 147615,
 'Aequorea victoria': 6100,
 'Aequorea victoria ': 6100,
 'Aequorea victoria  ': 6100,
 'Aequoria victoria': 6100,
 'Agaricia agaricites': 89882,
 'Agaricia fragilis': 165097,
 'Anemonia majano': 105399,
 'Anemonia rustica': 444856,
 'Anemonia sulcata': 6108,
 'Anemonia sulcata ': 6108,
 'Anthomedusae': 406427,
 'Anthomedusae sp.': None,
 'Anthomedusae sp. ': None,
 'Anthomedusae sp. DC-2005': 328397,
 'Astrangia lajollaensis': 262533,
 'Branchiostoma floridae': 7739,
 'Branchiostoma floridae ': 7739,
 'Branchiostoma lanceolatum': 7740,
 'Catalaphyllia jardinei': 46748,
 'Ceriantharia': 37512,
 'Cerianthus membranaceus': 208460,
 'Cerianthus sp.': 51771,
 'Chiridius poppei': 286301,
 'Clavularia sp.': 86521,
 'Clavularia sp. ': 86521,
 'Clytia gregaria': 27801,
 'Clytia hemisphaerica': 252671,
 'Cnidopus japonicus': 58804,
 'Condylactis gigantea': 47073,
 'Condylactis gigantea ': 47073,
 'Condylactis passiflora ': 175772,
 'Corynactis californica': 44298,
 'Cyphastrea microphthalma': 570133,
 'Danafungia horrida': 486202,
 'Dendronephthya sp.': 51110,
 'Dendronephthya sp. ': 51110,
 'Discosoma': 86599,
 'Discosoma sp': 86600,
 'Discosoma sp.': 86600,
 'Discosoma sp. ': 86600,
 'Discosoma sp. .': 86600,
 'Discosoma sp. 3': 86600,
 'Discosoma sp. LW-2004': 301246,
 'Discosoma sp.1': 86600,
 'Discosoma sp.2': 86600,
 'Discosoma striata': 105400,
 'Echinophyllia  ': 126651,
 'Echinophyllia echinata': 351729,
 'Echinophyllia sp.': None,
 'Echinophyllia sp. SC22': 301887,
 'Echinopora Sp. ': None,
 'Echinopora forskaliana': 526284,
 'Echinopora sp. ': None,
 'Entacmaea quadricolor': 6118,
 'Eusmilia fastigiata': 214972,
 'Favia favus': 102203,
 'Favia favus ': 102203,
 'Favites abdita': 126655,
 'Favites complanata': 498483,
 'Fungia concinna': 496660,
 'Galaxea fascicularis': 46745,
 'Geniopora djiboutiensis': 351727,
 'Goniopora tenuidens': 75301,
 'Herbaspirillum frisingense': 92645,
 'Heteractis Crispa': 175771,
 'Heteractis crispa': 175771,
 'Heteractis crispa ': 175771,
 'Heteractis magnifica ': 38281,
 'Hydnophora rigida': 46740,
 'Labidocera aestiva ': 163467,
 'Lobactis scutaria': 46714,
 'Lobophyllia hemprichii': 46758,
 'Lobophyllia hemprichii ': 46758,
 'Lobophyllia hemprichii  ': 46758,
 'Meandrina meandrites': 51056,
 'Merulina sp. ': None,
 'Monastraea cavernosa': 63558,
 'Montastraea annularis': 48500,
 'Montastraea cavernosa': 63558,
 'Montastraea faveolata': 48498,
 'Montastrea cavernosa': 63558,
 'Montipora efflorescens': 105610,
 'Montipora efflorescens ': 105610,
 'Montipora millepora': 351731,
 'Montipora sp.': None,
 'Montipora sp. ': None,
 'Montipora sp.  ': None,
 'Mycedium elephantotus': 51060,
 'Obelia sp.': 70918,
 'Pectiniidae': 46733,
 'Pectiniidae ': 46733,
 'Pectiniidae sp.  ': None,
 'Phialidium sp.': 1955689,
 'Phialidium sp. ': 1955689,
 'Platygira lamellina ': 242771,
 'Pocillopora damicornis': 46731,
 'Pontella meadi': 239965,
 'Pontella meadi ': 239965,
 'Pontella mimocerami': 661578,
 'Pontellidae sp.': None,
 'Pontellina plumata ': 239963,
 'Porites porites': 104760,
 'Psammocora sp.': None,
 'Psammocora superficialis': 371657,
 'Ptilosarcus sp.': None,
 'Renilla Reniformis': 6136,
 'Renilla muelleri': 37510,
 'Renilla reniformis': 6136,
 'Ricordea florida': 165100,
 'Sarcophyton sp.': None,
 'Sarcophyton sp. ': None,
 'Scleractinia sp. ': 1913369,
 'Scolymia cubensis': 165099,
 'Sphingomonas sp.': 28214,
 'Stylocoeniella sp. ': None,
 'Stylophora pistillata': 50429,
 'Symphyllia sp.': None,
 'Synthetic Construct': 32630,
 'Trachyphyllia geoffroyi': 196280,
 'Verrillofungia concinna': 496660,
 'Verrillofungia concinna  ': 496660,
 'Verrillofungia concinna  (Fungia concinna)': 496660,
 'Verrillofungia concinna (Fungia concinna)': 496660,
 'Zoanthus sp.': 105402,
 'Zoanthus sp. ': 105402,
 'Zoanthus sp.2': 105402}

BASEDIR = os.path.dirname(os.path.dirname(__file__))

def name_check(name):
    switch = None
    if 'before and after' in name.lower():
        if 'activation' in name:
            switch = 'pa'
        elif 'conversion' in name:
            switch = 'pc'
        elif 'switching' in name:
            switch = 'ps'
        else:
            switch = 'o'
        name = name.split(' (')[0]
    elif 'oxidized and reduced' in name.lower():
        switch = True
    return name.strip(), switch


def parensplit(st):
    return st.strip(')').replace('(', '').split()


def import_fpd(file=None):
    if file and os.path.isfile(os.path.join(BASEDIR, '_data', file)):
        file = os.path.join(BASEDIR, '_data', file)
    if not file:
        file = os.path.join(BASEDIR, '_data/FPD.csv')
    data = tablib.Dataset().load(open(file).read())

    errors = []
    for rownum, row in enumerate(data.dict):
        try:
            print(row['name'])
            # add the doi to the database
            # this is a workaround to use the web form, which asks for "reference_doi"
            # but converts it to the primary_reference primary key in the form view
            # try:
            #     doi = row.get('reference_doi', '').lstrip('https://dx.doi.org/')
            #     if doi:
            #         r, c = Reference.objects.get_or_create(doi=doi)
            #     if r:
            #         row['reference_doi'] = r.doi
            #     if c:
            #         print("created reference {}".format(r))
            # except Exception:
            #     r = None
            # add organism to database

            # just remove bad sequences
            try:
                protein_sequence_validator(row.get('seq'))
            except Exception:
                row['seq'] = None
            row['agg'] = row['agg'].strip()

            org = None
            if ORGLOOKUP.get(row.get('parent_organism')):
                org, ocreated = Organism.objects.get_or_create(id=ORGLOOKUP[row['parent_organism']])
                if ocreated:
                    print("created organism {}".format(org))
            row['parent_organism'] = org.pk if org else None

            # look for photoswitching in naming
            row['name'], switch = name_check(row['name'])

            namemismatch = False  # will be true if something already has this sequence with a diff name

            # check if protein already exists a variety of ways
            if row.get('seq') and Protein.objects.filter(seq=row.get('seq')).count():
                p = Protein.objects.get(seq=row.get('seq'))
            elif Protein.objects.filter(slug=slugify(row.get('name'))).count():
                p = Protein.objects.get(slug=slugify(row.get('name')))
            elif Protein.objects.filter(genbank=row['genbank']).count():
                p = Protein.objects.get(genbank=row['genbank'])
            else:
                p = None

            # don't overwrite existing values hack...
            if p:
                if p.genbank:
                    if row['genbank']:
                        if not p.genbank == row['genbank']:
                            pass
                            # errors.append('GenBank mismatch between {} and {}'.format(p.name, row.get('name')))
                    else:
                        row['genbank'] = p.genbank
                if p.uniprot:
                    if row['uniprot']:
                        if not p.uniprot == row['uniprot']:
                            pass
                            # errors.append('UniProt mismatch between {} and {}'.format(p.name, row.get('name')))
                    else:
                        row['uniprot'] = p.uniprot
                if not p.name == row.get('name', None):
                    namemismatch = True
                    # errors.append('same sequence name mismatch between {} and {}'.format(p.name, row.get('name')))
                    row['aliases'] = row.get('name')
                    row['name'] = p.name

            # create the protein form and validate/sve
            pform = forms.ProteinUpdateForm(row, instance=p) if p else forms.ProteinSubmitForm(row)
            pform.fields['reference_doi'].required = False
            if pform.is_valid():
                p = pform.save()
                doi = pform.cleaned_data.get('reference_doi')
                if pform.cleaned_data.get('reference_doi'):
                    ref, created = Reference.objects.get_or_create(doi=doi)
                    if created:
                        print("created Reference {}".format(ref))
                    ref.proteins.add(p)
                    if not p.primary_reference:
                        p.primary_reference = ref
                        p.save()
            else:
                errors.append("name: {}, row: {}, {}".format(data.dict[rownum]['name'], rownum, pform.errors.as_text()))

            # if we still don't have a protein instance, just move on
            if not p:
                print("STILL NO PROTEIN")
                continue

            if switch:
                states = [row.copy(), row.copy()]
                states[0]['name'] = 'before'
                states[1]['name'] = 'after'
                for col in ('ex_max', 'em_max', 'qy', 'ext_coeff', 'pka'):
                    if row[col]:
                        vals = parensplit(row[col])
                        for i, v in enumerate(vals):
                            try:
                                states[i][col] = float(v) if col in ('qy', 'pka') else int(v)
                            except Exception:
                                states[i][col] = None
                    else:
                        states[0][col] = None
                        states[1][col] = None

                if switch == 'pa' and not states[0]['em_max']:
                    transwave = states[0]['em_max']
                    states[0]['ex_max'] = None
                    states[0]['em_max'] = None
                    states[0]['is_dark'] = True
                else:
                    transwave = None
            else:
                states = [row.copy()]
                if not namemismatch:
                    states[0]['name'] = 'default'

            sinstances = []
            for snum, state in enumerate(states):
                if not state.get('is_dark', 0) and not (state.get('ex_max', 0) and state.get('em_max', 0)):
                    continue
                try:
                    if State.objects.filter(protein=p,
                            ex_max=state['ex_max'], em_max=state['em_max'],
                            ext_coeff=state['ext_coeff'], qy=state['qy']).count():
                        print('skipping already imported state on {}'.format(p.name))
                        continue
                except Exception:
                    pass
                state['protein'] = p.pk
                sform = forms.StateSubmitForm(state)
                if sform.is_valid():
                    sinstances.append(sform.save())
                else:
                    errors.append("name: {}, state: {}, {}".format(data.dict[rownum]['name'], snum, sform.errors.as_text()))

            if switch and len(sinstances) > 1:
                try:
                    StateTransition.objects.create(
                        protein=p,
                        from_state=sinstances[0],
                        to_state=sinstances[1],
                        transwave=transwave if transwave else None,
                    )
                    if switch == 'ps':
                        StateTransition.objects.create(
                            protein=p,
                            from_state=sinstances[1],
                            to_state=sinstances[0],
                        )
                except Exception:
                    print('failed to link states {} > {} '.format())

            p.save()  # register states

        except KeyboardInterrupt:
            raise
        except Exception as ex:
            print(ex)

    return errors


@require_superuser
def reload_all(seqs=False):
    importCSV()
    importPSFPs()
    importSeqs()
    importSpectra()
    import_organisms()
    if seqs:
        for P in Protein.objects.all():
            Q = fetch_ipg_sequence(P.name)
            if Q:
                P.seq = Q[1]
                P.ipg_id = Q[0]
                P.save()


#########################
#  Conversion
#########################

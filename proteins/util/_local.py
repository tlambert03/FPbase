import os
import json
import tablib
import traceback
import pandas as pd
from fpbase.users.models import User
from django.utils.text import slugify
import re
from references.models import Reference
from .. import forms
from ..models import (Protein, State, StateTransition, BleachMeasurement,
                      Organism, ProteinCollection, Dye, Spectrum)
from ..forms import SpectrumForm
from ..util.importers import text_to_spectra, import_chroma_spectra
from ..util.helpers import zip_wave_data
from ..util.spectra_import import import_spectral_data
from ..validators import protein_sequence_validator


try:
    SUPERUSER = User.objects.filter(is_superuser=True).first()
except Exception:
    SUPERUSER = None


BASEDIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


def get_nonan(obj, item):
    val = obj.get(item, None)
    if val is not None:
        if pd.isna(val):
            return None
    return val


def add_ref_to_prot(protein, doi, showexisting=False):
    ref, created = Reference.objects.get_or_create(doi=doi.strip())
    rf = 1 if created else 0
    protein.references.add(ref)
    if not protein.primary_reference:
        protein.primary_reference = ref
    protein.save()
    return rf


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


def importPSFPs(file=None):
    if file is None:
        url = 'https://raw.githubusercontent.com/kthorn/FPvisualization/master/PSFPs.csv'
        df = pd.read_csv(url)
    else:
        df = pd.read_csv(file)

    aggLookup = {
        'Tetramer': 't',
        'Monomer': 'm',
        'Dimer': 'd',
        'Dimer/Trimer': 'd',
    }

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
                    doi = re.sub('^https?://(dx\.)?doi.org/', '', prot.Source)
                    rf += add_ref_to_prot(p, doi)
            except Exception as e:
                # traceback.print_exc()
                print('Error importing reference: {}'.format(e))

    print("{} Sequences added; {} References imported".format(sq, rf))


def create_collection(name='FPvis Collection', desc='Proteins selected by Kurt Thorn at fpvis.org'):
    url = 'https://raw.githubusercontent.com/FPvisualization/fpvisualization.github.io/master/FPs.csv'
    df = pd.read_csv(url)
    col = ProteinCollection.objects.create(name=name, description=desc, owner=SUPERUSER)
    col.save()
    for n in df['Name']:
        try:
            p = Protein.objects.get(name=n)
            col.proteins.add(p)
        except Exception:
            print('{} failed'.format(n))
            pass


def import_csv_spectra(file, **kwargs):
    ''' import CSV or text file of spectral data

    kwargs:
    headers=None
    categories=[]  {d, p, l, f, c}
    stypes=[]  {ex, ab, em, 2p, bp, bx, bm, sp, lp, bs, qe, pd}
    owner=None
    minmax=None

    '''
    if not os.path.isfile(file):
        raise FileNotFoundError('Cannot find file: {}'.format(file))
    with open(file, 'r') as f:
        text = f.read()
    waves, data, headers = text_to_spectra(text)
    return import_spectral_data(waves, data, headers, **kwargs)


def import_thermo():
    d = os.path.join(BASEDIR, '_data/Thermo')
    for f in os.listdir(d):
        if not f.endswith('.csv'):
            continue
        name = f.strip('.csv')
        objs, errs = import_csv_spectra(os.path.join(d, f), categories=Spectrum.DYE, owner=name)
        if len(objs):
            owner = objs[0].owner
            for spect in objs:
                if spect.subtype == Spectrum.EX:
                    owner.ex_max = spect.peak_wave
                elif spect.subtype == Spectrum.EM:
                    owner.em_max = spect.peak_wave
                elif spect.subtype == Spectrum.TWOP:
                    owner.twop_ex_max = spect.peak_wave
            owner.manufacturer = 'ThermoFisher'
            part = name.lower().replace(' ', '-')
            owner.part = part
            owner.created_by = User.objects.first()
            owner.save()
        elif len(errs):
            print('Error importing: ', f)
            print(errs[0][1].as_text())

    update_dyes()


def update_dyes(file=None):
    if not file:
        file = os.path.join(BASEDIR, '_data/dyes.csv')
    if os.path.isfile(file):
        import tablib
        with open(os.path.join(BASEDIR, '_data/dyes.csv')) as csvf:
            D = tablib.Dataset().load(csvf.read())
            for row in D.dict:
                try:
                    owner = Dye.objects.get(slug=row['slug'])
                    for field in ('ext_coeff', 'qy', 'lifetime', 'pka', 'part', 'url'):
                        if row[field]:
                            try:
                                setattr(owner, field, float(row[field]))
                            except ValueError:
                                setattr(owner, field, row[field])
                    owner.save()
                except Exception as e:
                    print('Skipping {}: {}'.format(row['slug'], e))


def dyes_csv(file='/Users/talley/Desktop/dyes.csv'):
    import csv
    with open(file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # write your header first
        headers = ('id', 'slug', 'name', 'part', 'manufacturer', 'qy', 'ext_coeff', 'pka', 'lifetime')
        writer.writerow(headers)
        for obj in Dye.objects.all():
            row = []
            for field in headers:
                row.append(getattr(obj, field))
            writer.writerow(row)


def import2P():
    d = os.path.join(BASEDIR, '_data/2p_spectra')

    for f in os.listdir(d):
        name = f.replace('_w.txt', '')
        qs = Protein.objects.filter(name__iexact=name)
        if qs.exists():
            P = qs.first()
            infile = os.path.join(d, f)

            if not os.path.isfile(infile):
                raise FileNotFoundError('Cannot find file: {}'.format(infile))
            with open(infile, 'r') as f:
                text = f.read()

            x, y, headers = text_to_spectra(text)
            D = zip_wave_data(x, y[0])

            sf = SpectrumForm({
                'data': D,
                'category': Spectrum.PROTEIN,
                'subtype': Spectrum.TWOP,
                'owner': P.name
            })
            if sf.is_valid():
                obj = sf.save()
                P.default_state.twop_ex_max = obj._peakwave2p
                P.default_state.twop_peakGM = obj._peakval2p
                P.default_state.save()

                # add drobizhev reference
                ref, created = Reference.objects.get_or_create(doi='10.1038/nmeth.1596')
                P.references.add(ref)
                P.save()
                print('Successfuly import 2P spectrum for {}'.format(P.name))
            else:
                print('error on {}'.format(P.name))
                print(sf.errors.as_text())


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


def import_fpd(file=None, overwrite=True):
    if file and os.path.isfile(os.path.join(BASEDIR, '_data', file)):
        file = os.path.join(BASEDIR, '_data', file)
    if not file:
        file = os.path.join(BASEDIR, '_data/FPD.csv')
    data = tablib.Dataset().load(open(file).read())

    errors = []
    for rownum, row in enumerate(data.dict):
        try:
            for k, v in row.items():
                row[k] = v.strip()
            print(row['name'])

            # just remove bad sequences
            try:
                protein_sequence_validator(row.get('seq'))
            except Exception:
                row['seq'] = None
            row['agg'] = row['agg'].strip()

            if row.get('pdb', False):
                row['pdb'] = [i.strip() for i in row['pdb'].split(',') if i.strip()]
            else:
                row['pdb'] = []

            if row.get('aliases', False):
                row['aliases'] = [i.strip() for i in row['aliases'].split(',') if i.strip()]
            else:
                row['aliases'] = []

            org = None
            if row.get('parent_organism'):
                org, ocreated = Organism.objects.get_or_create(id=row['parent_organism'])
                if ocreated:
                    print("created organism {}".format(org))
            row['parent_organism'] = org.pk if org else None

            # look for photoswitching in naming
            row['name'], switch = name_check(row['name'])

            namemismatch = False  # will be true if something already has this sequence with a diff name

            # check if protein already exists a variety of ways
            if Protein.objects.filter(slug=slugify(row.get('name'))).exists():
                p = Protein.objects.get(slug=slugify(row.get('name')))
                if Protein.objects.filter(seq=row.get('seq')).exists():
                    p2 = Protein.objects.get(seq=row.get('seq'))
                    if p2 != p:
                        row['seq'] = None
                        errors.append('cannot assign {} sequence, {} already has it'.format(row.get('name'), p2.name))
            elif row.get('seq') and Protein.objects.filter(seq=row.get('seq')).exists():
                p = Protein.objects.get(seq=row.get('seq'))
            elif Protein.objects.filter(genbank=row['genbank']).exists():
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
                    print('Protein "{}" already has the same sequence as {}... skipping'.format(p.name, row['name']))
                    continue
                    # namemismatch = True
                    # errors.append('same sequence name mismatch between {} and {}'.format(p.name, row.get('name')))
                    # row['aliases'].append(row.get('name'))
                    # row['name'] = p.name

            # create the protein form and validate/sve
            ref = None
            pform = forms.ProteinForm(row, instance=p) if p else forms.ProteinForm(row)
            if pform.is_valid():
                p = pform.save() if overwrite else pform.save_new_only()
                doi = pform.cleaned_data.get('reference_doi')
                if pform.cleaned_data.get('reference_doi'):
                    ref, created = Reference.objects.get_or_create(doi=doi)
                    if created:
                        print("created Reference {}".format(ref))
                    ref.proteins.add(p)
                    if not p.primary_reference:
                        p.primary_reference = ref
                        p = pform.save() if overwrite else pform.save_new_only()
                if row['additional_refs']:
                    for doi in row['additional_refs'].split(','):
                        ref, created = Reference.objects.get_or_create(doi=doi.strip())
                        ref.proteins.add(p)
                        if created:
                            print("created Reference {}".format(ref))
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
                sform = forms.StateForm(state)
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

            if row['bleach']:
                state = sinstances[-1] if len(sinstances) else None
                bm = BleachMeasurement(rate=float(row['bleach']), reference=ref, state=state)
                bm.save()

            p = pform.save() if overwrite else pform.save_new_only()  # register states

        except KeyboardInterrupt:
            raise
        except Exception as ex:
            print(ex)

    return errors


def importMutations():
    from proteins.validators import validate_mutation

    file = os.path.join(BASEDIR, '_data/merged.csv')

    data = tablib.Dataset().load(open(file).read())
    mutOut = tablib.Dataset()
    mutOut.headers = ['name', 'parent', 'mutations', 'seq', 'doi']

    for n, row in enumerate(data.dict):
        if row['mutation']:
            mut = [m.strip() for m in row['mutation'].split('/') if m]
            name = row['name']
            parent = ''
            try:
                for i, m in enumerate(mut):
                    try:
                        validate_mutation(m)
                    except Exception:
                        if i == 0:
                            parent = mut.pop(0)
                        else:
                            raise
            except Exception:
                print(i, n)
                print("Failed:              ", name)
            mutOut.append((name, parent, "/".join(mut), row['seq'], row['reference_doi']))


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
 'Zoanthus sp.2': 105402,
 'Anthomedusae sp.': 406427,
 'Acanthastrea sp.': 406427,
 }


def import_chroma():
    parts = ('ZET488/594m', 'ZET488/561m', 'ZET405/488/561/647x', 'ZET405/488/561/647m',
             'ZET405/488/561/640m', 'T660lpxr', 'T570lp', 'T515lp', 'T455lp', 'T425lpxr',
             'T400lp', 'S630/60m', 'S535/40m', 'S470/30m', 'Q660lp', 'Q585lp', 'Q565lp',
             'Q505lp', 'HQ630/40m', 'HQ535/50m', 'HQ480/40x', 'ET705/72m', 'ET700/75m',
             'ET645/30x', 'ET632/60m', 'ET620/60x', 'ET620/60m', 'ET610lp',
             'ET605/70m', 'ET605/52m', 'ET600/50m', 'ET572/35x', 'ET560/40x', 'ET555/25x',
             'ET545/30x', 'ET535/50m', 'ET535/30m', 'ET525/50m', 'ET525/36m',
             'ET510/80m', 'ET500lp', 'ET500/20x', 'ET490/20x', 'ET480/40x',
             'ET480/40m', 'ET480/30x', 'ET470/24m', 'ET460/50m',
             'ET455/50m', 'ET436/20x', 'ET430/24x', 'ET402/15x', 'ET395/25x', 'ET380x',
             'ET340x', 'D620/60m', 'D540/25x', 'D/F/Cy3/Cy5', 'CFP/YFP/mCherry XT',
             'AT350/50x', '89100bs', '86002v1bs', '69008bs', '69002bs', '59022bs')

    for p in parts:
        try:
            import_chroma_spectra(p)
        except Exception as e:
            print('Could not import chroma part {} ({})'.format(p, e))


def import_lights():
    file = os.path.join(BASEDIR, '_data/broadband_light_spectra.csv')
    objs, errs = import_csv_spectra(file,
                                    categories=Spectrum.LIGHT,
                                    stypes=Spectrum.PD)
    for obj in objs:
        obj.owner.created_by = User.objects.first()
        obj.created_by = User.objects.first()
        obj.save()
        obj.owner.save()


def import_lumencor():
    D = os.path.join(BASEDIR, '_data/lumencor')
    for f in os.listdir(D):
        if not f.endswith('.txt'):
            continue
        objs, errs = import_csv_spectra(os.path.join(D, f),
                                        categories=Spectrum.LIGHT,
                                        stypes=Spectrum.PD,
                                        owner=f.replace('.txt', ''))
        if errs:
            print(errs[0][1].as_text())
        if objs:
            owner = objs[0].owner
            owner.manufacturer = 'lumencor'
            owner.part = slugify(f.strip('.txt'))
            if 'spectrax' in f.lower():
                owner.url = 'http://lumencor.com/products/spectra-x-light-engine/'
            owner.created_by = User.objects.first()
            owner.save()

    D = os.path.join(BASEDIR, '_data/lumencorFilters')
    for f in os.listdir(D):
        if not f.endswith('.txt'):
            continue
        name = f.strip('.txt')
        name = name[:3] + '/' + name[-2:] + 'x'
        ownername = 'Lumencor ' + name
        objs, errs = import_csv_spectra(os.path.join(D, f),
                                        categories=Spectrum.FILTER,
                                        stypes=Spectrum.BPX,
                                        owner=ownername)
        if errs:
            print(errs[0][1].as_text())
        if objs:
            owner = objs[0].owner
            owner.manufacturer = 'lumencor'
            owner.part = slugify(name)
            owner.url = 'http://lumencor.com/products/filters-for-spectra-x-light-engines/'
            owner.created_by = User.objects.first()
            owner.save()

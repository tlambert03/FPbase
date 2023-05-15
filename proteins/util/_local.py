import csv
import json
import os
import re
import traceback

import pandas as pd
import requests
import tablib
from django.utils.text import slugify

from fpbase.users.models import User
from references.models import Reference

from .. import forms
from ..forms import SpectrumForm
from ..models import (
    BleachMeasurement,
    Dye,
    Organism,
    OSERMeasurement,
    Protein,
    ProteinCollection,
    Spectrum,
    State,
    StateTransition,
)
from ..util.helpers import getprot, zip_wave_data
from ..util.importers import import_chroma_spectra, text_to_spectra
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
    """
    mainly intended as an import function for the proteins in the NIC table
    """
    if file is None:
        url = "https://raw.githubusercontent.com/FPvisualization/fpvisualization.github.io/master/FPs.csv"
        df = pd.read_csv(url)
    else:
        df = pd.read_csv(file)

    ps = 0
    st = 0
    rf = 0
    for _i, prot in df.iterrows():
        if not Protein.objects.filter(slug=slugify(prot.Name)).exists():
            print(f"importing {prot.Name}...")
            p = Protein(
                name=prot.Name,
                created_by=SUPERUSER,
                agg=get_nonan(prot, "agg"),
                switch_type=Protein.BASIC,
            )
            p.save()
            ps += 1

            s = State(
                # name        = 'default',  # now the default state name
                ex_max=get_nonan(prot, "lambda_ex"),
                em_max=get_nonan(prot, "lambda_em"),
                ext_coeff=get_nonan(prot, "E"),
                qy=get_nonan(prot, "QY"),
                pka=get_nonan(prot, "pka"),
                maturation=get_nonan(prot, "mature"),
                lifetime=get_nonan(prot, "lifetime"),
                protein=p,
                created_by=SUPERUSER,
            )
            s.save()
            st += 1

            # add it as default state
            p.default_state = s
            p.save()

            # add bleach numbers
            if get_nonan(prot, "bleach"):
                BleachMeasurement.objects.update_or_create(rate=get_nonan(prot, "bleach"), state=s)
        else:
            print(f"{prot.Name} already in database...")
            p = Protein.objects.get(slug=slugify(prot.Name))

        try:
            if pd.isna(prot.DOI) or not prot.DOI.startswith("10"):
                continue
            for doi in prot.DOI.split(" "):
                if not (doi and doi.startswith("10")):
                    continue
                rf += add_ref_to_prot(p, doi)
        except Exception as e:
            # traceback.print_exc()
            print(f"error importing reference: {e}")

    print(f"{ps} Proteins, {st} States, and {rf} References imported")


def linkstates(df):
    linksurl = "https://raw.githubusercontent.com/FPvisualization/fpvisualization.github.io/master/links.csv"
    linksdf = pd.read_csv(linksurl)
    q = df.set_index("UID").to_dict()

    for _i, link in linksdf.iterrows():
        p = Protein.objects.get(slug=slugify(q["Name"][link.state1]))
        fromState = p.states.get(name=q["state"][link.state1])
        toState = p.states.get(name=q["state"][link.state2])

        t, created = StateTransition.objects.get_or_create(
            protein=p,
            trans_wave=int(link.lambda_sw),
            from_state=fromState,
            to_state=toState,
        )
        if created:
            print(f"created: {t}")


def importPSFPs(file=None):
    if file is None:
        url = "https://raw.githubusercontent.com/kthorn/FPvisualization/master/PSFPs.csv"
        df = pd.read_csv(url)
    else:
        df = pd.read_csv(file)

    aggLookup = {"Tetramer": "t", "Monomer": "m", "Dimer": "d", "Dimer/Trimer": "d"}

    ps = 0
    st = 0
    rf = 0
    for _i, prot in df.iterrows():
        try:
            p, created = Protein.objects.get_or_create(
                slug=slugify(prot.Name),
                defaults={
                    "name": prot.Name,
                    "created_by": SUPERUSER,
                    "switch_type": prot.type,
                    "agg": aggLookup[prot.Aggregation] if get_nonan(prot, "Aggregation") else None,
                },
            )
            if created:
                print(f"PROTEIN CREATED: {prot.Name}")
                ps += 1
            else:
                print(f"Protein Found  : {prot.Name}")

            if not pd.isna(prot.DOI) and prot.DOI.startswith("10"):
                try:
                    rf += add_ref_to_prot(p, prot.DOI)
                except Exception as e:
                    # traceback.print_exc()
                    print(f"Error importing reference: {e}")

            state, created = p.states.get_or_create(
                name=prot.state,
                defaults={
                    "name": prot.state,
                    "is_dark": "True" if "off" in prot.state.lower() else False,
                    "ex_max": get_nonan(prot, "lambda_ex"),
                    "em_max": get_nonan(prot, "lambda_em"),
                    "ext_coeff": get_nonan(prot, "E"),
                    "qy": get_nonan(prot, "QY"),
                    "pka": get_nonan(prot, "pka"),
                    "maturation": get_nonan(prot, "mature"),
                    "lifetime": get_nonan(prot, "lifetime"),
                    "protein": p,
                    "created_by": SUPERUSER,
                },
            )

            # needs work... would like ON state to show up in tables...
            # but only if it is a DARK -> ON protein
            if p.switch_type == "ps" or p.switch_type == "pa":
                if not state.is_dark:
                    p.default_state = state
                    p.save()
            elif prot.get("initialState", False):
                p.default_state = state
                p.save()

            if get_nonan(prot, "bleach"):
                BleachMeasurement.objects.update_or_create(rate=get_nonan(prot, "bleach"), state=state)
            if created:
                print(f"STATE CREATED  : {prot.state}")
                st += 1
            else:
                print(f"State Found    : {prot.state}")

        except Exception as e:
            traceback.print_exc()
            print(f"failed to import {prot.Name}: {e}")

    print(f"{ps} Proteins, {st} States, and {rf} References imported")
    linkstates(df)


def importSeqs(file=None):
    if file is None:
        basedir = os.path.dirname(os.path.dirname(__file__))
        url = os.path.join(basedir, "_data/FPseqs.csv")
        df = pd.read_csv(url)
    else:
        df = pd.read_csv(file)

    rf = 0
    sq = 0
    for _i, prot in df.iterrows():
        if Protein.objects.filter(name__icontains=prot.Name).count() == 1:
            p = Protein.objects.get(name__icontains=prot.Name)
            if p.seq is None:
                p.seq = prot.AminoAcidSequence
                p.save()
                print(f"Added sequence to {prot.Name}")
                sq += 1
            else:
                seq = prot.AminoAcidSequence.upper()
                seq = "".join(seq.split())
                if p.seq.upper() != seq:
                    print(f"Non-matching sequence found for {prot.Name}!")
            try:
                if "dx.doi" in prot.Source:
                    doi = re.sub(r"^https?://(dx\.)?doi.org/", "", prot.Source)
                    rf += add_ref_to_prot(p, doi)
            except Exception as e:
                # traceback.print_exc()
                print(f"Error importing reference: {e}")

    print(f"{sq} Sequences added; {rf} References imported")


def create_collection(name="FPvis Collection", desc="Proteins selected by Kurt Thorn at fpvis.org"):
    url = "https://raw.githubusercontent.com/FPvisualization/fpvisualization.github.io/master/FPs.csv"
    df = pd.read_csv(url)
    col = ProteinCollection.objects.create(name=name, description=desc, owner=SUPERUSER)
    col.save()
    for n in df["Name"]:
        try:
            p = Protein.objects.get(name=n)
            col.proteins.add(p)
        except Exception:
            print(f"{n} failed")
            pass


def import_csv_spectra(file, **kwargs):
    """import CSV or text file of spectral data

    kwargs:
    headers=None
    categories=[]  {d, p, l, f, c}
    stypes=[]  {ex, ab, em, 2p, bp, bx, bm, sp, lp, bs, qe, pd}
    owner=None
    minmax=None

    """
    if not os.path.isfile(file):
        raise FileNotFoundError(f"Cannot find file: {file}")
    with open(file) as f:
        text = f.read()
    waves, data, headers = text_to_spectra(text)
    return import_spectral_data(waves, data, headers, **kwargs)


def import_thermo():
    d = os.path.join(BASEDIR, "_data/Thermo")
    for f in os.listdir(d):
        if not f.endswith(".csv"):
            continue
        name = f.strip(".csv")
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
            owner.manufacturer = "ThermoFisher"
            part = name.lower().replace(" ", "-")
            owner.part = part
            owner.created_by = User.objects.first()
            owner.save()
        elif len(errs):
            print("Error importing: ", f)
            print(errs[0][1].as_text())

    update_dyes()


def import_atto():
    from django.core.files import File

    from proteins.forms import SpectrumForm

    d = os.path.join(BASEDIR, "_data/ATTO")
    props = pd.DataFrame.from_csv(os.path.join(BASEDIR, "_data/ATTO/_attoprops.csv"))
    props = props.T.to_dict()

    for f in os.listdir(d):
        if not f.endswith(".txt"):
            continue
        try:
            name = f.split("_")[0].replace("ATTO", "ATTO ")
            with open(os.path.join(d, f), "rb") as _f:
                djfile = File(_f)
                sf = SpectrumForm(
                    {
                        "owner": name,
                        "category": Spectrum.DYE,
                        "subtype": Spectrum.ABS if "abs.txt" in f else Spectrum.EM,
                        "ph": 7.4,
                        "solven": "PBS",
                    },
                    {"file": djfile},
                )
                if sf.is_valid():
                    obj = sf.save(commit=False)
                    obj.created_by = User.objects.first()
                    obj.save()

                    dye = obj.owner
                    dye.manufacturer = "ATTO-TEC GmbH"
                    dye.part = name
                    uri = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/{}.pdf"
                    url = uri.format(name.replace(" ", "_"))
                    if requests.get(url).status_code == 200:
                        dye.url = url
                    if name in props:
                        for k, v in props[name].items():
                            if v:
                                setattr(dye, k, v)
                    else:
                        print("could not find props for ", name)
                    dye.created_by = User.objects.first()
                    dye.save()

                    print("SUCCESS: ", obj)
                else:
                    print(sf.errors)
        except Exception:
            raise


def import_oser(file=None):
    if not file:
        file = os.path.join(BASEDIR, "_data/oser.csv")
    if os.path.isfile(file):
        import tablib

        with open(file) as csvf:
            D = tablib.Dataset().load(csvf.read())
            for row in D.dict:
                try:
                    prot = Protein.objects.get(name=row["name"])
                    ref = Reference.objects.get(doi=row["doi"])
                    OSERMeasurement.objects.create(
                        protein=prot,
                        reference=ref,
                        percent=float(row.get("percent")) if row.get("percent", False) else None,
                        percent_stddev=float(row.get("percent_stddev")) if row.get("percent_stddev", False) else None,
                        percent_ncells=int(row.get("percent_ncells")) if row.get("percent_ncells", False) else None,
                        oserne=float(row.get("oserne")) if row.get("oserne", False) else None,
                        oserne_stddev=float(row.get("oserne_stddev")) if row.get("oserne_stddev", False) else None,
                        oserne_ncells=int(row.get("oserne_ncells")) if row.get("oserne_ncells", False) else None,
                        celltype=row.get("celltype", None),
                    )
                except Exception as e:
                    print("Skipping {}: {}".format(row["name"], e))


def update_dyes(file=None):
    if not file:
        file = os.path.join(BASEDIR, "_data/dyes.csv")
    if os.path.isfile(file):
        import tablib

        with open(os.path.join(BASEDIR, "_data/dyes.csv")) as csvf:
            D = tablib.Dataset().load(csvf.read())
            for row in D.dict:
                try:
                    owner = Dye.objects.get(slug=row["slug"])
                    for field in ("ext_coeff", "qy", "lifetime", "pka", "part", "url"):
                        if row[field]:
                            try:
                                setattr(owner, field, float(row[field]))
                            except ValueError:
                                setattr(owner, field, row[field])
                    owner.save()
                except Exception as e:
                    print("Skipping {}: {}".format(row["slug"], e))


def dyes_csv(file="/Users/talley/Desktop/dyes.csv"):
    import csv

    with open(file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # write your header first
        headers = (
            "id",
            "slug",
            "name",
            "part",
            "manufacturer",
            "qy",
            "ext_coeff",
            "pka",
            "lifetime",
        )
        writer.writerow(headers)
        for obj in Dye.objects.all():
            row = []
            for field in headers:
                row.append(getattr(obj, field))
            writer.writerow(row)


def import2P():
    d = os.path.join(BASEDIR, "_data/2p_spectra")

    for f in os.listdir(d):
        name = f.replace("_w.txt", "")
        qs = Protein.objects.filter(name__iexact=name)
        if qs.exists():
            P = qs.first()
            infile = os.path.join(d, f)

            if not os.path.isfile(infile):
                raise FileNotFoundError(f"Cannot find file: {infile}")
            with open(infile) as f:
                text = f.read()

            x, y, headers = text_to_spectra(text)
            D = zip_wave_data(x, y[0])

            sf = SpectrumForm(
                {
                    "data": D,
                    "category": Spectrum.PROTEIN,
                    "subtype": Spectrum.TWOP,
                    "owner": P.name,
                }
            )
            if sf.is_valid():
                obj = sf.save()
                P.default_state.twop_ex_max = obj._peakwave2p
                P.default_state.twop_peakGM = obj._peakval2p
                P.default_state.save()

                # add drobizhev reference
                ref, created = Reference.objects.get_or_create(doi="10.1038/nmeth.1596")
                P.references.add(ref)
                P.save()
                print(f"Successfuly import 2P spectrum for {P.name}")
            else:
                print(f"error on {P.name}")
                print(sf.errors.as_text())


def name_check(name):
    switch = None
    if "before and after" in name.lower():
        if "activation" in name:
            switch = "pa"
        elif "conversion" in name:
            switch = "pc"
        elif "switching" in name:
            switch = "ps"
        else:
            switch = "o"
        name = name.split(" (")[0]
    elif "oxidized and reduced" in name.lower():
        switch = True
    return name.strip(), switch


def parensplit(st):
    return st.strip(")").replace("(", "").split()


def import_fpd(file=None, overwrite=True):
    if file and os.path.isfile(os.path.join(BASEDIR, "_data", file)):
        file = os.path.join(BASEDIR, "_data", file)
    if not file:
        file = os.path.join(BASEDIR, "_data/FPD.csv")
    data = tablib.Dataset().load(open(file).read())

    errors = []
    for rownum, row in enumerate(data.dict):
        try:
            for k, v in row.items():
                row[k] = v.strip()
            print(row["name"])

            # just remove bad sequences
            try:
                protein_sequence_validator(row.get("seq"))
            except Exception:
                row["seq"] = None
            row["agg"] = row["agg"].strip()

            if row.get("pdb", False):
                row["pdb"] = [i.strip() for i in row["pdb"].split(",") if i.strip()]
            else:
                row["pdb"] = []

            if row.get("aliases", False):
                row["aliases"] = [i.strip() for i in row["aliases"].split(",") if i.strip()]
            else:
                row["aliases"] = []

            org = None
            if row.get("parent_organism"):
                org, ocreated = Organism.objects.get_or_create(id=row["parent_organism"])
                if ocreated:
                    print(f"created organism {org}")
            row["parent_organism"] = org.pk if org else None

            # look for photoswitching in naming
            row["name"], switch = name_check(row["name"])

            namemismatch = False  # will be true if something already has this sequence with a diff name

            # check if protein already exists a variety of ways
            if Protein.objects.filter(slug=slugify(row.get("name"))).exists():
                p = Protein.objects.get(slug=slugify(row.get("name")))
                if Protein.objects.filter(seq=row.get("seq")).exists():
                    p2 = Protein.objects.get(seq=row.get("seq"))
                    if p2 != p:
                        row["seq"] = None
                        errors.append("cannot assign {} sequence, {} already has it".format(row.get("name"), p2.name))
            elif row.get("seq") and Protein.objects.filter(seq=row.get("seq")).exists():
                p = Protein.objects.get(seq=row.get("seq"))
            elif Protein.objects.filter(genbank=row["genbank"]).exists():
                p = Protein.objects.get(genbank=row["genbank"])
            else:
                p = None

            # don't overwrite existing values hack...
            if p:
                if p.genbank:
                    if row["genbank"]:
                        if not p.genbank == row["genbank"]:
                            pass
                            # errors.append('GenBank mismatch between {} and {}'.format(p.name, row.get('name')))
                    else:
                        row["genbank"] = p.genbank
                if p.uniprot:
                    if row["uniprot"]:
                        if not p.uniprot == row["uniprot"]:
                            pass
                            # errors.append('UniProt mismatch between {} and {}'.format(p.name, row.get('name')))
                    else:
                        row["uniprot"] = p.uniprot

                if not p.name == row.get("name", None):
                    print('Protein "{}" already has the same sequence as {}... skipping'.format(p.name, row["name"]))
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
                doi = pform.cleaned_data.get("reference_doi")
                if pform.cleaned_data.get("reference_doi"):
                    ref, created = Reference.objects.get_or_create(doi=doi)
                    if created:
                        print(f"created Reference {ref}")
                    ref.proteins.add(p)
                    if not p.primary_reference:
                        p.primary_reference = ref
                        p = pform.save() if overwrite else pform.save_new_only()
                if row["additional_refs"]:
                    for doi in row["additional_refs"].split(","):
                        ref, created = Reference.objects.get_or_create(doi=doi.strip())
                        ref.proteins.add(p)
                        if created:
                            print(f"created Reference {ref}")
            else:
                errors.append(
                    "name: {}, row: {}, {}".format(data.dict[rownum]["name"], rownum, pform.errors.as_text())
                )

            # if we still don't have a protein instance, just move on
            if not p:
                print("STILL NO PROTEIN")
                continue

            if switch:
                states = [row.copy(), row.copy()]
                states[0]["name"] = "before"
                states[1]["name"] = "after"
                for col in ("ex_max", "em_max", "qy", "ext_coeff", "pka"):
                    if row[col]:
                        vals = parensplit(row[col])
                        for i, v in enumerate(vals):
                            try:
                                states[i][col] = float(v) if col in ("qy", "pka") else int(v)
                            except Exception:
                                states[i][col] = None
                    else:
                        states[0][col] = None
                        states[1][col] = None

                if switch == "pa" and not states[0]["em_max"]:
                    transwave = states[0]["em_max"]
                    states[0]["ex_max"] = None
                    states[0]["em_max"] = None
                    states[0]["is_dark"] = True
                else:
                    transwave = None
            else:
                states = [row.copy()]
                if not namemismatch:
                    states[0]["name"] = "default"

            sinstances = []
            for snum, state in enumerate(states):
                if not state.get("is_dark", 0) and not (state.get("ex_max", 0) and state.get("em_max", 0)):
                    continue
                try:
                    if State.objects.filter(
                        protein=p,
                        ex_max=state["ex_max"],
                        em_max=state["em_max"],
                        ext_coeff=state["ext_coeff"],
                        qy=state["qy"],
                    ).count():
                        print(f"skipping already imported state on {p.name}")
                        continue
                except Exception:
                    pass
                state["protein"] = p.pk
                sform = forms.StateForm(state)
                if sform.is_valid():
                    sinstances.append(sform.save())
                else:
                    errors.append(
                        "name: {}, state: {}, {}".format(data.dict[rownum]["name"], snum, sform.errors.as_text())
                    )

            if switch and len(sinstances) > 1:
                try:
                    StateTransition.objects.create(
                        protein=p,
                        from_state=sinstances[0],
                        to_state=sinstances[1],
                        transwave=transwave if transwave else None,
                    )
                    if switch == "ps":
                        StateTransition.objects.create(protein=p, from_state=sinstances[1], to_state=sinstances[0])
                except Exception:
                    print("failed to link states")

            if row["bleach"]:
                state = sinstances[-1] if len(sinstances) else None
                bm = BleachMeasurement(rate=float(row["bleach"]), reference=ref, state=state)
                bm.save()

            p = pform.save() if overwrite else pform.save_new_only()  # register states

        except KeyboardInterrupt:
            raise
        except Exception as ex:
            print(ex)

    return errors


def importMutations():
    from proteins.validators import validate_mutation

    file = os.path.join(BASEDIR, "_data/merged.csv")

    data = tablib.Dataset().load(open(file).read())
    mutOut = tablib.Dataset()
    mutOut.headers = ["name", "parent", "mutations", "seq", "doi"]

    for n, row in enumerate(data.dict):
        if row["mutation"]:
            mut = [m.strip() for m in row["mutation"].split("/") if m]
            name = row["name"]
            parent = ""
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
            mutOut.append((name, parent, "/".join(mut), row["seq"], row["reference_doi"]))


def import_organisms():
    with open("_data/species.json") as f:
        D = json.load(f)
    for k in D.keys():
        o = Organism(id=k)
        o.save()


ORGLOOKUP = {
    "Acanthastrea sp. ": None,
    "Acropora aculeus": 287157,
    "Acropora digitifera": 70779,
    "Acropora eurostoma": 526283,
    "Acropora hyacinthus": 55974,
    "Acropora millepora": 45264,
    "Acropora nobilis": 70781,
    "Acropora pulchra": 140239,
    "Acropora sp.": None,
    "Acropora tenuis": 70783,
    "Actinia equina": 6106,
    "Actinia equina ": 6106,
    "Aequorea Victoria": 6100,
    "Aequorea Victoria ": 6100,
    "Aequorea Victoria  ": 6100,
    "Aequorea coerulescens": 210840,
    "Aequorea coerulescens ": 210840,
    "Aequorea macrodactyla": 147615,
    "Aequorea victoria": 6100,
    "Aequorea victoria ": 6100,
    "Aequorea victoria  ": 6100,
    "Aequoria victoria": 6100,
    "Agaricia agaricites": 89882,
    "Agaricia fragilis": 165097,
    "Anemonia majano": 105399,
    "Anemonia rustica": 444856,
    "Anemonia sulcata": 6108,
    "Anemonia sulcata ": 6108,
    "Anthomedusae": 406427,
    "Anthomedusae sp. ": None,
    "Anthomedusae sp. DC-2005": 328397,
    "Astrangia lajollaensis": 262533,
    "Branchiostoma floridae": 7739,
    "Branchiostoma floridae ": 7739,
    "Branchiostoma lanceolatum": 7740,
    "Catalaphyllia jardinei": 46748,
    "Ceriantharia": 37512,
    "Cerianthus membranaceus": 208460,
    "Cerianthus sp.": 51771,
    "Chiridius poppei": 286301,
    "Clavularia sp.": 86521,
    "Clavularia sp. ": 86521,
    "Clytia gregaria": 27801,
    "Clytia hemisphaerica": 252671,
    "Cnidopus japonicus": 58804,
    "Condylactis gigantea": 47073,
    "Condylactis gigantea ": 47073,
    "Condylactis passiflora ": 175772,
    "Corynactis californica": 44298,
    "Cyphastrea microphthalma": 570133,
    "Danafungia horrida": 486202,
    "Dendronephthya sp.": 51110,
    "Dendronephthya sp. ": 51110,
    "Discosoma": 86599,
    "Discosoma sp": 86600,
    "Discosoma sp.": 86600,
    "Discosoma sp. ": 86600,
    "Discosoma sp. .": 86600,
    "Discosoma sp. 3": 86600,
    "Discosoma sp. LW-2004": 301246,
    "Discosoma sp.1": 86600,
    "Discosoma sp.2": 86600,
    "Discosoma striata": 105400,
    "Echinophyllia  ": 126651,
    "Echinophyllia echinata": 351729,
    "Echinophyllia sp.": None,
    "Echinophyllia sp. SC22": 301887,
    "Echinopora Sp. ": None,
    "Echinopora forskaliana": 526284,
    "Echinopora sp. ": None,
    "Entacmaea quadricolor": 6118,
    "Eusmilia fastigiata": 214972,
    "Favia favus": 102203,
    "Favia favus ": 102203,
    "Favites abdita": 126655,
    "Favites complanata": 498483,
    "Fungia concinna": 496660,
    "Galaxea fascicularis": 46745,
    "Geniopora djiboutiensis": 351727,
    "Goniopora tenuidens": 75301,
    "Herbaspirillum frisingense": 92645,
    "Heteractis Crispa": 175771,
    "Heteractis crispa": 175771,
    "Heteractis crispa ": 175771,
    "Heteractis magnifica ": 38281,
    "Hydnophora rigida": 46740,
    "Labidocera aestiva ": 163467,
    "Lobactis scutaria": 46714,
    "Lobophyllia hemprichii": 46758,
    "Lobophyllia hemprichii ": 46758,
    "Lobophyllia hemprichii  ": 46758,
    "Meandrina meandrites": 51056,
    "Merulina sp. ": None,
    "Monastraea cavernosa": 63558,
    "Montastraea annularis": 48500,
    "Montastraea cavernosa": 63558,
    "Montastraea faveolata": 48498,
    "Montastrea cavernosa": 63558,
    "Montipora efflorescens": 105610,
    "Montipora efflorescens ": 105610,
    "Montipora millepora": 351731,
    "Montipora sp.": None,
    "Montipora sp. ": None,
    "Montipora sp.  ": None,
    "Mycedium elephantotus": 51060,
    "Obelia sp.": 70918,
    "Pectiniidae": 46733,
    "Pectiniidae ": 46733,
    "Pectiniidae sp.  ": None,
    "Phialidium sp.": 1955689,
    "Phialidium sp. ": 1955689,
    "Platygira lamellina ": 242771,
    "Pocillopora damicornis": 46731,
    "Pontella meadi": 239965,
    "Pontella meadi ": 239965,
    "Pontella mimocerami": 661578,
    "Pontellidae sp.": None,
    "Pontellina plumata ": 239963,
    "Porites porites": 104760,
    "Psammocora sp.": None,
    "Psammocora superficialis": 371657,
    "Ptilosarcus sp.": None,
    "Renilla Reniformis": 6136,
    "Renilla muelleri": 37510,
    "Renilla reniformis": 6136,
    "Ricordea florida": 165100,
    "Sarcophyton sp.": None,
    "Sarcophyton sp. ": None,
    "Scleractinia sp. ": 1913369,
    "Scolymia cubensis": 165099,
    "Sphingomonas sp.": 28214,
    "Stylocoeniella sp. ": None,
    "Stylophora pistillata": 50429,
    "Symphyllia sp.": None,
    "Synthetic Construct": 32630,
    "Trachyphyllia geoffroyi": 196280,
    "Verrillofungia concinna": 496660,
    "Verrillofungia concinna  ": 496660,
    "Verrillofungia concinna  (Fungia concinna)": 496660,
    "Verrillofungia concinna (Fungia concinna)": 496660,
    "Zoanthus sp.": 105402,
    "Zoanthus sp. ": 105402,
    "Zoanthus sp.2": 105402,
    "Anthomedusae sp.": 406427,
    "Acanthastrea sp.": 406427,
}


def import_chroma():
    parts = (
        "ZET488/594m",
        "ZET488/561m",
        "ZET405/488/561/647x",
        "ZET405/488/561/647m",
        "ZET405/488/561/640m",
        "T660lpxr",
        "T570lp",
        "T515lp",
        "T455lp",
        "T425lpxr",
        "T400lp",
        "S630/60m",
        "S535/40m",
        "S470/30m",
        "Q660lp",
        "Q585lp",
        "Q565lp",
        "Q505lp",
        "HQ630/40m",
        "HQ535/50m",
        "HQ480/40x",
        "ET705/72m",
        "ET700/75m",
        "ET645/30x",
        "ET632/60m",
        "ET620/60x",
        "ET620/60m",
        "ET610lp",
        "ET605/70m",
        "ET605/52m",
        "ET600/50m",
        "ET572/35x",
        "ET560/40x",
        "ET555/25x",
        "ET545/30x",
        "ET535/50m",
        "ET535/30m",
        "ET525/50m",
        "ET525/36m",
        "ET510/80m",
        "ET500lp",
        "ET500/20x",
        "ET490/20x",
        "ET480/40x",
        "ET480/40m",
        "ET480/30x",
        "ET470/24m",
        "ET460/50m",
        "ET455/50m",
        "ET436/20x",
        "ET430/24x",
        "ET402/15x",
        "ET395/25x",
        "ET380x",
        "ET340x",
        "D620/60m",
        "D540/25x",
        "D/F/Cy3/Cy5",
        "CFP/YFP/mCherry XT",
        "AT350/50x",
        "89100bs",
        "86002v1bs",
        "69008bs",
        "69002bs",
        "59022bs",
    )

    for p in parts:
        try:
            import_chroma_spectra(p)
        except Exception as e:
            print(f"Could not import chroma part {p} ({e})")


def import_lights():
    file = os.path.join(BASEDIR, "_data/broadband_light_spectra.csv")
    objs, errs = import_csv_spectra(file, categories=Spectrum.LIGHT, stypes=Spectrum.PD)
    for obj in objs:
        obj.owner.created_by = User.objects.first()
        obj.created_by = User.objects.first()
        obj.save()
        obj.owner.save()


def import_lumencor():
    D = os.path.join(BASEDIR, "_data/lumencor")
    for f in os.listdir(D):
        if not f.endswith(".txt"):
            continue
        objs, errs = import_csv_spectra(
            os.path.join(D, f),
            categories=Spectrum.LIGHT,
            stypes=Spectrum.PD,
            owner=f.replace(".txt", ""),
        )
        if errs:
            print(errs[0][1].as_text())
        if objs:
            owner = objs[0].owner
            owner.manufacturer = "lumencor"
            if f.endswith(".txt"):
                f = f[-4:]
            owner.part = slugify(f)
            if "spectrax" in f.lower():
                owner.url = "http://lumencor.com/products/spectra-x-light-engine/"
            owner.created_by = User.objects.first()
            owner.save()

    D = os.path.join(BASEDIR, "_data/lumencorFilters")
    for f in os.listdir(D):
        if not f.endswith(".txt"):
            continue
        name = f[-4:]
        name = name[:3] + "/" + name[-2:] + "x"
        ownername = "Lumencor " + name
        objs, errs = import_csv_spectra(
            os.path.join(D, f),
            categories=Spectrum.FILTER,
            stypes=Spectrum.BPX,
            owner=ownername,
        )
        if errs:
            print(errs[0][1].as_text())
        if objs:
            owner = objs[0].owner
            owner.manufacturer = "lumencor"
            owner.part = slugify(name)
            owner.url = "http://lumencor.com/products/filters-for-spectra-x-light-engines/"
            owner.created_by = User.objects.first()
            owner.save()


def osfp_import():
    from collections import defaultdict

    with open("_data/osfp-full-data-set.csv") as f:
        csvrows = csv.reader(f)
        D = defaultdict(dict)
        for name, agg, seq, doi in csvrows:
            D[name]["seq"] = seq.replace("\n", "")
            D[name]["agg"] = agg
            D[name]["doi"] = doi
        return D
        """
        for name, agg, seq, doi in csvrows:
            seq = seq.replace('\n', '')
            p = getname(name)
            if not p:
                continue
            if not p.seq:
                print('ADD SEQ: ', name)
            elif p.seq == seq:
                pass
            else:
                print('seq mismatch: ', name)

            # DOI
            if p.primary_reference:
                if not doi == p.primary_reference.doi:
                    print('DOI mismatch on {} ({}->{})'.format(name, p.primary_reference.doi, doi))
            elif doi:
                print('Add DOI to {}: {}'.format(p, doi))

            # AGG
            if not agg == p.get_agg_display():
                if not p.agg == Protein.WEAK_DIMER:
                    print('change {} agg from {} to {}'.format(name, p.get_agg_display(), agg))
        """


def snapgene_import():
    with open("snapgene.csv") as f:
        csvrows = csv.reader(f)
        for row in csvrows:
            name = row[4]
            # seq = row[6]
            try:
                p = Protein.objects.get(name=name)
            except Protein.DoesNotExist:
                continue
            if not p:
                continue
            if p.seq:
                pass
                # if not ParasailAlignment.from_seqs(seq, p.seq).mutations:
                #    continue
                # print(name, ParasailAlignment.from_seqs(seq, p.seq).mutations)


def get_gb_data(file):
    with open(file) as handle:
        text = handle.read()
    q = ("DEFINITION", "ACCESSION", "VERSION", "KEYWORDS", "SOURCE")
    pat = ""
    for i in range(len(q) - 1):
        pat += r"(?=.*{} (?P<{}>.+){})?".format(q[i], q[i].lower().replace(" ", ""), q[i + 1])
    pat += r"(?=.*COMMENT (?P<comment>.+)FEAT)?"
    pat += r"(?=.*PUBMED\s+(?P<pub>.+)REF)?"
    pat += r'(?=.*/translation="(?P<seq>.+)")?'
    D = re.search(pat, text, re.DOTALL).groupdict()
    for k, v in D.items():
        if v:
            D[k] = v.replace("\n", "").replace("  ", "").strip()
    D["seq"] = D.get("seq", "").replace(" ", "")
    D["name"] = os.path.basename(file).strip(".gb")
    return D


def import_tree(filepath=None):
    from proteins.models import Lineage, Protein, Reference
    from proteins.validators import validate_mutationset

    if not filepath:
        filepath = os.path.join(BASEDIR, "_data/Lineage.xlsx")

    data = pd.read_excel(filepath)
    data = data.where((pd.notnull(data)), "").astype(str)

    nonames = []
    for (
        _rownum,
        (doi, _author, _year, name, parnt, _mutation, aliases, _note),
    ) in data.iterrows():
        if name:
            try:
                getprot(name)
            except Protein.DoesNotExist:
                nonames.append((name, doi, parnt, aliases))
    if nonames:
        # r = input('The following {} names do not exist in the database:\n{}'
        #          '\nWould you like to create them? (y/n):  '
        #          .format(len(nonames), "\n".join([x[0] for x in nonames])))
        # if r.lower() == 'y':
        while True:
            i = 0
            n = 0
            for name, doi, parnt, aliases in nonames:
                try:
                    if parnt:
                        par = getprot(parnt)
                        org = par.parent_organism
                    else:
                        org = None
                except Exception:
                    org = None
                if doi:
                    ref, _ = Reference.objects.get_or_create(doi=doi.lower())
                else:
                    ref = None

                if not Protein.objects.filter(name=name).exists():
                    p = Protein.objects.create(
                        name=name,
                        created_by=SUPERUSER,
                        primary_reference=ref,
                        parent_organism=org,
                    )
                    print("Created new protein: ", p)
                    if aliases:
                        p.aliases = [a.strip() for a in aliases.split(",")]
                    p.status = "approved"
                    p.save()
                    n += 1

            i += 1
            if i > 100:
                print("more than 100 iterations...")
                break
            if n == 0:
                break

    for doi in data["ref"]:
        if doi:
            doi = doi.lower()
            try:
                Reference.objects.get(doi=doi)
            except Exception:
                try:
                    print("creating reference doi: ", doi)
                    Reference.objects.create(doi__iexact=doi, created_by=SUPERUSER)
                except Exception as e:
                    print(f"failed to add reference {doi}: {e}")

    print("starting ... ")
    count = 1
    while count > 0:
        count = 0
        for (
            _rownum,
            (doi, _author, _year, prot, parnt, mutation, _alias, _note),
        ) in data.iterrows():
            prot = prot.strip()
            validate_mutationset(mutation)
            if not len(prot):
                continue
            try:
                child = getprot(prot)
            except Protein.DoesNotExist:
                print(f'Could not find child "{prot}"')
                continue
            parent = None
            try:
                parnt = parnt.strip()
                parent = getprot(parnt)
            except Protein.DoesNotExist:
                if parnt:
                    print("Could not find parent ", parnt)
                    continue
            try:
                ref = Reference.objects.get(doi__iexact=doi)
            except Reference.DoesNotExist:
                ref = None

            if parent and not Lineage.objects.filter(protein=child).exists():
                try:
                    parent = Lineage.objects.get(protein=parent)
                except Lineage.DoesNotExist:
                    continue
                L = Lineage.objects.create(
                    protein=child,
                    parent=parent,
                    reference=ref,
                    mutation=mutation,
                    created_by=SUPERUSER,
                )
                print(f"Created {L.parent} -> {L.protein}")
                count += 1
            elif child and not Lineage.objects.filter(protein=child).exists():
                # create a root node
                L = Lineage.objects.create(
                    protein=child,
                    parent=None,
                    reference=ref,
                    mutation=mutation,
                    created_by=SUPERUSER,
                )

                print(f"Created root: {L.protein}")
                count += 1
            # except IntegrityError as e:
            #     # already created this one...
            #     # print("\t\tIntegrityError: {}".format(e))
            #     pass
            # except Lineage.DoesNotExist as e:
            #     # parent doesn't exist yet... wait a minute
            #     print("DoesNotExist: {}".format(e))
            #     pass
            # except Exception as e:
            #     print(e)
            #     try:
            #         print('\tParent: ', parent)
            #     except Exception:
            #         pass
            #     print('\tChild: ', child)
            #     print(e)
            #     raise

    # add missing parent orgs:
    for name, _, _, _ in nonames:
        if Lineage.objects.filter(protein__name=name).exists:
            L = Lineage.objects.get(protein__name=name)
            root = L.get_root()
            p = getprot(name)
            p.parent_organism = root.protein.parent_organism
            if p.name.startswith("d") and not p.agg:
                p.agg = "d"
            if p.name.startswith("m") and not p.agg:
                p.agg = "m"
            if p.name.startswith("td") and not p.agg:
                p.agg = "td"
            p.save()
    from proteins.util.maintain import add_missing_seqs

    add_missing_seqs()
    add_missing_seqs()
    add_missing_seqs()
    [s.save() for s in Lineage.objects.all()]


CF_EC = {
    "CF350": 18000,
    "CF405S": 33000,
    "CF405M": 41000,
    "CF405L": 24000,
    "CF430": 40000,
    "CF440": 40000,
    "CF450": 40000,
    "CF488A": 70000,
    "CF514": 105000,
    "CF532": 96000,
    "CF535ST": 95000,
    "CF543": 100000,
    "CF555": 150000,
    "CF568": 100000,
    "CF594": 115000,
    "CF594ST": 115000,
    "CF620R": 115000,
    "CF633": 100000,
    "CF640R": 105000,
    "CF647": 240000,
    "CF660C": 100000,
    "CF680": 210000,
    "CF680R": 140000,
    "CF750": 250000,
    "CF770": 220000,
    "CF790": 210000,
    "CF800": 210000,
}


def import_cf(path="_data/Biotium/biotium_spectra.csv"):
    from proteins.forms import SpectrumForm

    df = pd.DataFrame.from_csv(os.path.join(BASEDIR, path))

    for dye_name in df:
        data = zip_wave_data(df.index, df[dye_name])
        subtype = Spectrum.ABS if dye_name.endswith("Abs") else Spectrum.EM
        name = dye_name.replace(" Abs", "").replace(" Em", "")
        sf = SpectrumForm(
            {
                "owner": name,
                "category": Spectrum.DYE,
                "subtype": subtype,
                "data": data,
                "confirmation": True,
            }
        )
        if sf.is_valid():
            obj = sf.save(commit=False)
            obj.created_by = User.objects.first()
            obj.save()

            dye = obj.owner
            dye.manufacturer = "Biotium"
            dye.part = name
            # dye.url = 'https://biotium.com/'
            uri = (
                "https://biotium.com/search-results/keyword/{}/search-in/product/cat-in"
                + "/all/search-other/product,p_sku,p_cat,post,page"
            )
            dye.url = uri.format(name.replace(" ", "+"))
            # if name in props:
            #     for k, v in props[name].items():
            #         if v:
            #             setattr(dye, k, v)
            # else:
            #     print('could not find props for ', name)
            if name in CF_EC:
                dye.ext_coeff = CF_EC[name]
            dye.created_by = User.objects.first()
            dye.save()

            print("SUCCESS: ", obj)
        else:
            print(sf.errors)

import csv
import os
import re

import pandas as pd
from django.utils.text import slugify

from fpbase.users.models import User

from ..models import (
    Protein,
    Spectrum,
)
from ..util.helpers import getprot, zip_wave_data
from ..util.importers import import_chroma_spectra, text_to_spectra
from ..util.spectra_import import import_spectral_data

try:
    SUPERUSER = User.objects.filter(is_superuser=True).first()
except Exception:
    SUPERUSER = None


BASEDIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


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
    objs, _errs = import_csv_spectra(file, categories=Spectrum.LIGHT, stypes=Spectrum.PD)
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

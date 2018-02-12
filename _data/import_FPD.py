from tablib import Dataset
from proteins.validators import protein_sequence_validator

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


def get_diffs(rowa, rowb, headers):
    difflist = []
    for col, vala in enumerate(rowa):
        valb = rowb[col]
        if vala and valb and vala != valb:
            difflist.append(headers[col])
    return difflist


FPs_file = '/Users/talley/Dropbox (HMS)/Python/fpbase/_data/FPs.csv'
FPs = Dataset().load(open(FPs_file).read())
FPdict = {}
for row in FPs.dict:
    row['reference_doi'] = row['reference_doi'].strip()
    FPdict[row['name']] = row

file = '/Users/talley/Dropbox (HMS)/Python/fpbase/_data/FPD.csv'
data = Dataset().load(open(file).read())
parens = Dataset()
noparens = Dataset()
parens.headers = data.headers
noparens.headers = data.headers

# general data cleaning and sorting
ni = data.headers.index('name')
oi = data.headers.index('parent_organism')
di = data.headers.index('reference_doi')
si = data.headers.index('seq')
for row in data:
    r = [r.strip().replace('\n', ' ') for r in row]
    r[oi] = ORGLOOKUP.get(r[oi], None)
    r[di] = r[di].lstrip('http://dx.doi.org/').strip()
    try:
        protein_sequence_validator(r[si])
    except Exception:
        r[si] = None
    if r[ni].endswith(')'):
        parens.append(r)
    else:
        noparens.append(r)

# ## FIND unique rows and duplicated rows ###
noparens_uniq = Dataset()
noparens_uniq.headers = data.headers
noparens_dupe = Dataset()
noparens_dupe.headers = data.headers
for ind, row in enumerate(noparens):
    if noparens['name'].count(row[0]) > 1:
        noparens_dupe.append(row)
    else:
        noparens_uniq.append(row)


# clean up duplicates
namedupes = set()
for row in noparens_dupe:
    namedupes.add(tuple([i for i, x in enumerate(noparens_dupe['name']) if x == row[0]]))

dellist = []
for dupepair in namedupes:
    if len(dupepair) == 2:
        rowA, rowB = noparens_dupe[dupepair[0]], noparens_dupe[dupepair[1]]
        diffs = get_diffs(rowA, rowB, data.headers)
        if not any([i in diffs for i in ('seq', 'ex_max', 'em_max')]):
            newrow = [a or b for a, b in zip(rowA, rowB)]
            for d in diffs:
                if d == 'reference_doi':
                    id = data.headers.index(d)
                    ia = data.headers.index('additional_refs')
                    if rowA[ia]:
                        newrow[ia] = rowA[ia] + ', ' + rowB[id]
                    else:
                        newrow[ia] = rowB[id]
                else:
                    i = data.headers.index(d)
                    newrow[i] = rowA[i] + ', ' + rowB[i]
            noparens_uniq.append(newrow)
            dellist.append(rowA[0])

noparens_dupe_temp = Dataset()
noparens_dupe_temp.headers = data.headers
for row in noparens_dupe:
    if not row[0] in dellist:
        noparens_dupe_temp.append(row)
noparens_dupe = noparens_dupe_temp


# pull out the multi-state proteins
switchers = Dataset()
switchers.headers = data.headers
D = []
for i, row in enumerate(parens):
    phrases = ('Before and After', 'pH', '(Red', '(Green', 'State', 'Form', 'IrisFP')
    if any([phrase in row[0] for phrase in phrases]):
        switchers.append(row)
        D.append(i)
D.reverse()
for i in D:
    del parens[i]


# with open('/Users/talley/Desktop/fps.csv', 'w') as file:
#    file.write(FPs.export('csv'))
with open('/Users/talley/Desktop/parens.csv', 'w') as file:
    file.write(parens.export('csv'))
with open('/Users/talley/Desktop/switchers.csv', 'w') as file:
    file.write(switchers.export('csv'))
# with open('/Users/talley/Desktop/noparens_uniq.csv', 'w') as file:
#    file.write(noparens_uniq.export('csv'))
with open('/Users/talley/Desktop/noparens_dupe.csv', 'w') as file:
    file.write(noparens_dupe.export('csv'))


# look for conflicting info
for row in noparens_uniq.dict:
    exclude = []
    n = row['name']
    if n in FPdict:
        e = []
        for k in data.headers:
            if k in exclude:
                continue
            if row[k] and FPdict[n][k] and not FPdict[n][k] == row[k]:
                e.append(k)
        if len(e):
            print("{}: {}".format(n, e))


# merge good stuff
merged = Dataset().load(open(FPs_file).read())
for row in noparens_uniq:
    n = row[0]
    if n in merged['name']:
        popped = merged[merged['name'].index(n)]
        del merged[merged['name'].index(n)]
        comborow = tuple([a or b for a, b in zip(popped, row)])
        merged.append(comborow)
    else:
        if not row[di]:  # skip without reference
            continue
        merged.append(row)

merged = merged.sort(0)


#oser = Dataset().load(open('/Users/talley/Dropbox (HMS)/Python/fpbase/_data/oser.csv').read())



with open('/Users/talley/Desktop/merged.csv', 'w') as file:
    file.write(merged.export('csv'))


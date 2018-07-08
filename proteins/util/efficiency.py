import numpy as np


def spectral_product(arrlist):
    # calculate (overlapping) product of a list of spectra.values
    # assumes monotonic increase with step = 1

    minwaves = [int(a[0][0]) for a in arrlist]
    minwave = max(minwaves)
    offsets = [minwave - i for i in minwaves]
    waverange = min([int(a[-1][0]) for a in arrlist]) - minwave
    vals = []
    for w in range(waverange + 1):
        vals.append([w + minwave, np.product([a[w + offsets[i]][1]
                     for i, a in enumerate(arrlist)])])
    return vals


def area(arr):
    area = 0
    for i in range(len(arr) - 1):
        area += (arr[i][1] + arr[i + 1][1]) / 2
    return area


def print_scope_report(report, sortby='bright', show=('bright', 'ex', 'em'), limit=10):
    form = '{:<20}' + ('{:>15}  ' * limit)
    for oc, values in report.items():
        topfluors = sorted(values.items(), key=lambda kv: kv[1][sortby])
        topfluors.reverse()
        topfluors = topfluors[:limit]
        print(form.format(oc.name, *[t[0].fluor_name for t in topfluors]))
        for item in show:
            print(form.format(item, *[t[1][item] for t in topfluors]))
        print()


def microscope_efficiency_report(scope, *args, **kwargs):
    return oclist_efficiency_report(scope.optical_configs.all(), *args, **kwargs)


def oclist_efficiency_report(oclist, owner_collection):
    D = {}
    for oc in oclist:
        print("OPTICAL CONFIG: ", oc)
        D[oc] = path_efficiency_report(oc, owner_collection)
    return D


def path_efficiency_report(oc, owner_collection):
    # where oc is an optical config and owner_collection is a list of a SpectrumOwner subclass
    D = {}
    oc_em = spectral_product(oc.em_spectra)
    oc_ex = spectral_product(oc.ex_spectra)
    for owner in owner_collection:
        print("\t", owner)
        D[owner] = {
            'ex': None,
            'em': None,
            'bright': 0
        }
        if owner.em_spectrum:
            combospectrum = spectral_product([oc_em, owner.em_spectrum.data])
            D[owner]['em'] = round(area(combospectrum) / area(owner.em_spectrum.data), 3)
        if owner.ex_spectrum:
            combospectrum = spectral_product([oc_ex, owner.ex_spectrum.data])
            D[owner]['ex'] = round(area(combospectrum) / area(oc_ex), 3)
        if D[owner]['em'] and D[owner]['ex'] and owner.ext_coeff and owner.qy:
            b = D[owner]['em'] * D[owner]['ex'] * owner.ext_coeff * owner.qy / 1000
            D[owner]['bright'] = round(b, 3)
    return D

import csv

import numpy as np
from django.http import HttpResponse
from scipy import interpolate
from scipy.signal import (
    argrelextrema,  # savgol_filter
    savgol_filter,
)


def is_monotonic(array):
    array = np.array(array)
    return np.all(array[1:] > array[:-1])


def make_monotonic(x, y):
    x, y = list(zip(*sorted(zip(x, y))))
    x, xind = np.unique(x, return_index=True)
    y = np.array(y)[xind]
    return x, y


def interp_linear(x, y, s=1, savgol=False):
    """Interpolate pair of vectors at integer increments between min(x) and max(x)"""
    if not is_monotonic(x):
        x, y = make_monotonic(x, y)
    xnew = range(int(np.ceil(min(x))), int(np.floor(max(x))))
    F = interpolate.interp1d(x, y)
    ynew = F(xnew)
    if savgol:
        ynew = savgol_filter(ynew, 9, 2)
    return xnew, ynew


def interp_univar(x, y, s=1, savgol=False):
    if not is_monotonic(x):
        x, y = make_monotonic(x, y)
    """Interpolate pair of vectors at integer increments between min(x) and max(x)"""
    xnew = range(int(np.ceil(min(x))), int(np.floor(max(x))))
    f = interpolate.InterpolatedUnivariateSpline(x, y)
    ynew = f(xnew)
    if savgol:
        ynew = norm2one(savgol_filter(ynew, 15, 2))
    return xnew, ynew


def norm2one(y):
    return [round(max(yy / max(y), 0), 4) for yy in y]


def step_size(lol):
    x, y = zip(*lol)
    s = set(np.subtract(x[1:], x[:-1]))
    if len(s) != 1:  # multiple step sizes
        return 0
    return s.pop()


def norm2P(y):
    """Normalize peak value of vector to one"""
    y = np.array(y)
    localmax = argrelextrema(y, np.greater, order=100)
    # can't be within first 10 points
    localmax = [i for i in localmax[0] if i > 10]
    if not localmax:
        return y, 0, 0
    maxind = localmax[np.argmax(y[localmax])]
    maxy = y[maxind]
    return [round(max(yy / maxy, 0), 4) for yy in y], maxy, maxind


def spectra2csv(spectralist, filename="fpbase_spectra.csv"):
    globalmin = int(min(sp.min_wave for sp in spectralist))
    globalmax = int(max(sp.max_wave for sp in spectralist))
    response = HttpResponse(content_type="text/csv")
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    writer = csv.writer(response)
    sp_headers = []
    wave_dicts = {}
    for sp in spectralist:
        sp_headers += [str(sp)]
        wave_dicts[str(sp)] = dict(sp.data)
    writer.writerow(["wavelength", *sp_headers])
    for wave in range(globalmin, globalmax + 1):
        row = [wave]
        row += [wave_dicts[sp].get(wave, "") for sp in sp_headers]
        writer.writerow(row)
    return response


def interp2int(x, y, s=1):
    """Interpolate pair of vectors at integer increments between min(x) and max(x)"""
    xnew = range(int(min(x)), int(max(x)))
    # ynew = cubic_interp1d(xnew, x, y)  # if no scipy
    tck = interpolate.splrep(x, y, s=s)
    ynew = interpolate.splev(xnew, tck, der=0)
    ynew = savgol_filter(ynew, 15, 2)
    return xnew, ynew


if __name__ == "__main__":
    import json
    import subprocess
    import sys

    from proteins.util.importers import text_to_spectra

    def file2spectra(file, dtype="", getcol=0):
        waves, outdata, headers = text_to_spectra(file)
        x = waves
        y = outdata[getcol]
        spectra = [list(i) for i in zip(x, y)]
        if dtype.lower() == "json":
            import json

            spectra = json.dumps(spectra)
        return spectra

    def set_clipboard_data(data):
        data = json.dumps(data, separators=(",", ":"))
        p = subprocess.Popen(["pbcopy"], stdin=subprocess.PIPE)
        p.stdin.write(bytes(data.strip('"'), "utf-8"))
        p.stdin.close()
        p.wait()

    infile = sys.argv[1]
    with open(infile) as f:
        waves, outdata, headers = text_to_spectra(f.read())
    out = [[float(n) for n in a] for a in zip(waves, outdata[0])]
    set_clipboard_data(out)
    print("data copied to clipboard")

from scipy import interpolate
from scipy.signal import savgol_filter
from django.http import HttpResponse
import csv


def spectra2csv(spectralist, filename='fpbase_spectra.csv'):
    globalmin = int(min([sp.min_wave for sp in spectralist]))
    globalmax = int(max([sp.max_wave for sp in spectralist]))
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)
    sp_headers = []
    wave_dicts = dict()
    for sp in spectralist:
        sp_headers += [str(sp)]
        wave_dicts[str(sp)] = dict(sp.data)
    writer.writerow(['wavelength'] + sp_headers)
    for wave in range(globalmin, globalmax + 1):
        row = [wave]
        row += [wave_dicts[sp].get(wave, '') for sp in sp_headers]
        writer.writerow(row)
    return response


def interp2int(x, y, s=1):
    '''Interpolate pair of vectors at integer increments between min(x) and max(x)'''
    xnew = range(int(min(x)), int(max(x)))
    # ynew = cubic_interp1d(xnew, x, y)  # if no scipy
    tck = interpolate.splrep(x, y, s=s)
    ynew = interpolate.splev(xnew, tck, der=0)
    ynew = savgol_filter(ynew, 15, 2)
    return xnew, ynew


if __name__ == '__main__':

    import json
    import subprocess
    import sys
    from proteins.util.importers import text_to_spectra

    def file2spectra(file, dtype='', getcol=0):
        waves, outdata, headers = text_to_spectra(file)
        x = waves
        y = outdata[getcol]
        spectra = [list(i) for i in zip(x, y)]
        if dtype.lower() == 'json':
            import json
            spectra = json.dumps(spectra)
        return spectra

    def setClipboardData(data):
        data = json.dumps(data, separators=(',', ':'))
        p = subprocess.Popen(['pbcopy'], stdin=subprocess.PIPE)
        p.stdin.write(bytes(data.strip('"'), 'utf-8'))
        p.stdin.close()
        p.wait()

    infile = sys.argv[1]
    with open(infile, 'r') as f:
        waves, outdata, headers = text_to_spectra(f.read())
    out = [list([float(n) for n in a]) for a in zip(waves, outdata[0])]
    setClipboardData(out)
    print('data copied to clipboard')


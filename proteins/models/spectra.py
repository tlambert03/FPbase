import tablib
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.signal import argrelextrema
import numpy as np

# if no scipy
# import numpy as np
# def cubic_interp1d(x0, x, y):
#     """
#     Interpolate a 1-D function using cubic splines.
#       x0 : a float or an 1d-array
#       x : (N,) array_like
#           A 1-D array of real/complex values.
#       y : (N,) array_like
#           A 1-D array of real values. The length of y along the
#           interpolation axis must be equal to the length of x.

#     Implement a trick to generate at first step the cholesky matrice L of
#     the tridiagonal matrice A (thus L is a bidiagonal matrice that
#     can be solved in two distinct loops).
#     """

#     x = np.asfarray(x)
#     y = np.asfarray(y)

#     # remove non finite values
#     # indexes = np.isfinite(x)
#     # x = x[indexes]
#     # y = y[indexes]

#     # check if sorted
#     if np.any(np.diff(x) < 0):
#         indexes = np.argsort(x)
#         x = x[indexes]
#         y = y[indexes]

#     size = len(x)

#     xdiff = np.diff(x)
#     ydiff = np.diff(y)

#     # allocate buffer matrices
#     Li = np.empty(size)
#     Li_1 = np.empty(size-1)
#     z = np.empty(size)

#     # fill diagonals Li and Li-1 and solve [L][y] = [B]
#     Li[0] = np.sqrt(2*xdiff[0])
#     Li_1[0] = 0.0
#     B0 = 0.0  # natural boundary
#     z[0] = B0 / Li[0]

#     for i in range(1, size-1, 1):
#         Li_1[i] = xdiff[i-1] / Li[i-1]
#         Li[i] = np.sqrt(2*(xdiff[i-1]+xdiff[i]) - Li_1[i-1] * Li_1[i-1])
#         Bi = 6*(ydiff[i]/xdiff[i] - ydiff[i-1]/xdiff[i-1])
#         z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

#     i = size - 1
#     Li_1[i-1] = xdiff[-1] / Li[i-1]
#     Li[i] = np.sqrt(2*xdiff[-1] - Li_1[i-1] * Li_1[i-1])
#     Bi = 0.0  # natural boundary
#     z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

#     # solve [L.T][x] = [y]
#     i = size-1
#     z[i] = z[i] / Li[i]
#     for i in range(size-2, -1, -1):
#         z[i] = (z[i] - Li_1[i-1]*z[i+1])/Li[i]

#     # find index
#     index = x.searchsorted(x0)
#     np.clip(index, 1, size-1, index)

#     xi1, xi0 = x[index], x[index-1]
#     yi1, yi0 = y[index], y[index-1]
#     zi1, zi0 = z[index], z[index-1]
#     hi1 = xi1 - xi0

#     # calculate cubic
#     f0 = zi0/(6*hi1)*(xi1-x0)**3 + \
#          zi1/(6*hi1)*(x0-xi0)**3 + \
#          (yi1/hi1 - zi1*hi1/6)*(x0-xi0) + \
#          (yi0/hi1 - zi0*hi1/6)*(xi1-x0)
#     return f0


def get_file_data(file):
    import re
    with open(file, 'r') as f:
        text = f.read()
        try:
            # find the first string, split on multiple tokens
            float(re.split(';|,|\n|\t', text)[0])
            data = tablib.Dataset().load(text, headers=[])
        except ValueError:
            # first number is not a float... probably a header
            data = tablib.Dataset().load(text)
    if data.width != 2:
        print("get_file_data only takes data from the first two columns")
    x = data['x'] if data.headers and 'x' in data.headers else data.get_col(0)
    y = data['y'] if data.headers and 'y' in data.headers else data.get_col(1)
    x = [float(xx) for xx in x]
    y = [float(yy) for yy in y]
    x, y = tuple(zip(*sorted(zip(x, y))))  # sort based on x
    return x, y


def interp2int(x, y, s=1):
    '''Interpolate pair of vectors at integer increments between min(x) and max(x)'''
    xnew = range(int(min(x)), int(max(x)))
    # ynew = cubic_interp1d(xnew, x, y)  # if no scipy
    tck = interpolate.splrep(x, y, s=s)
    ynew = interpolate.splev(xnew, tck, der=0)
    ynew = savgol_filter(ynew, 15, 2)
    return xnew, ynew


def interp_univar(x, y, s=1):
    '''Interpolate pair of vectors at integer increments between min(x) and max(x)'''
    xnew = range(int(min(x)), int(max(x)))
    F = interpolate.UnivariateSpline(x, y, s=s)
    ynew = F(xnew)
    # ynew = savgol_filter(ynew, 15, 2)
    return xnew, ynew


def interp_linear(x, y, s=1):
    '''Interpolate pair of vectors at integer increments between min(x) and max(x)'''
    xnew = range(int(min(x)), int(max(x)))
    F = interpolate.interp1d(x, y)
    ynew = F(xnew)
    ynew = savgol_filter(ynew, 9, 2)
    return xnew, ynew


def norm2one(y):
    '''Normalize peak value of vector to one'''
    return [round(max(yy/max(y), 0), 4) for yy in y]


def norm2P(y):
    '''Normalize peak value of vector to one'''
    localmax = argrelextrema(y, np.greater, order=100)
    # can't be within first 10 points
    localmax = [i for i in localmax[0] if i > 10]
    maxind = localmax[np.argmax(y[localmax])]
    maxy = y[maxind]
    return [round(max(yy/maxy, 0), 4) for yy in y], maxy, maxind


def file2spectra(file, interp=True, norm=True, dtype=''):
    x, y = get_file_data(file)
    if interp:
        x, y = interp2int(x, y)
    if norm:
        y = norm2one(y)
    spectra = [list(x) for x in zip(x, y)]
    if dtype.lower() == 'json':
        import json
        spectra = json.dumps(spectra)
    return spectra


# Run it as follows:
# x, y = get_file_data('/Users/talley/Desktop/CyOFP1em.txt')
# x, y = interp_and_norm(x, y)

if __name__ == '__main__':

    import json
    import subprocess
    import sys

    def setClipboardData(data):
        data = json.dumps(data, separators=(',', ':'))
        p = subprocess.Popen(['pbcopy'], stdin=subprocess.PIPE)
        p.stdin.write(bytes(data.strip('"'), 'utf-8'))
        p.stdin.close()
        p.wait()

    infile = sys.argv[1]
    x, y = get_file_data(infile)
    x, y = interp2int(x, y)
    y = norm2one(y)
    out = [list(a) for a in zip(x, y)]
    setClipboardData(out)
    print('data copied to clipboard')


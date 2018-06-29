import re
from uuid import uuid4


def shortuuid(padding=None):
    number = uuid4().int
    output = ""
    alph = list("23456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz")
    alpha_len = len(alph)
    while number:
        number, digit = divmod(number, alpha_len)
        output += alph[digit]
    if padding:
        remainder = max(padding - len(output), 0)
        output = output + alph[0] * remainder
    return output


def zip_wave_data(waves, data, minmax=None):
    minmax = minmax or (300, 1600)
    return [list(i) for i in zip(waves, data) if (i[1] > 0 and minmax[0] <= i[0] <= minmax[1])]


def wave_to_hex(wavelength, gamma=1):
    '''This converts a given wavelength into an approximate RGB value.
    The given wavelength is in nanometers.
    The range of wavelength is 380 nm through 750 nm.

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''
    if not wavelength:
        return '#000'

    wavelength = float(wavelength)
    if 520 <= wavelength:
        wavelength += 40

    if wavelength < 380:
        R = 0.05
        G = 0.0
        B = 0.15
    elif wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (770 - wavelength) / (770 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.18
        G = 0.0
        B = 0.05
    R *= 255
    G *= 255
    B *= 255
    return '#%02x%02x%02x' % (int(R), int(G), int(B))


#def wave_to_hex(wave):
#    wave = int(wave)
#    if wave < 380:
#        return "#000000"
#    if wave > 780:
#        return "#FFFF00"
#    else:
#        return COLORS[wave]

def get_color_group(ex_max, em_max):
    if (em_max - ex_max) > 90:
        return "Long Stokes Shift", "#80A0FF"
    if ex_max < 380:
        return "UV", "#C080FF"
    if ex_max < 421:
        return "Blue", "#8080FF"
    if ex_max < 473:
        return "Cyan", "#80FFFF"
    if ex_max < 505:
        return "Green", "#80FF80"
    if ex_max < 515:
        return "Green/Yellow", "#CCFF80"
    if ex_max < 531:
        return "Yellow", "#FFFF80"
    if ex_max < 555:
        return "Orange", "#FFC080"
    if ex_max < 600:
        return "Red", "#FFA080"
    if ex_max < 631:
        return "Far Red", "#FF8080"
    if ex_max < 800:
        return "Near IR", "#B09090"


def mless(name):
    if re.search('^m[A-Z]', name):
        return name.lstrip('m')
    if name.startswith('monomeric'):
        name = name.lstrip('monomeric')
    if name.startswith('Monomeric'):
        name = name.lstrip('Monomeric')
    return name.lstrip(' ')


def get_base_name(name):
    '''return core name of protein, stripping prefixes like "m" or "Tag"'''

    # remove PA/(Pa), PS, PC, from beginning
    if name.startswith(('PA', 'Pa', 'PS', 'Ps', 'PC', 'pc', 'rs')):
        name = name[2:]

    if re.match('LSS', name):
        name = name[3:].lstrip('-')

    # remove m (if next letter is caps) or monomeric
    if re.match('m[A-Z]', name):
        name = name[1:]

    # get rid of Td or td
    if re.match('[Tt][Dd][A-Z]', name):
        name = name[2:]

    name = name.lstrip('Monomeric')
    name = name.lstrip('Tag')

    # remove E at beginning (if second letter is caps)
    if re.match('E[A-Z]', name):
        name = name[1:]
    # remove S at beginning (if second letter is caps)
    if re.match('S[A-Z]', name):
        name = name[1:]
    # remove T- at beginning (if second letter is caps)
    if re.match('T-', name):
        name = name[2:]

    name = name.lstrip('-').lstrip(' ')

    return name


# ###########################################
#       Spectral Functions
# ###########################################

def calculate_spectral_overlap(donor, acceptor):
    accEx  = acceptor.default_state.ex_spectra
    accEC  = acceptor.default_state.ext_coeff
    donEm  = donor.default_state.em_spectra
    # donQY  = donor.default_state.qy
    donCum = sum(donEm.y)

    minAcc = accEx.min_wave
    maxAcc = accEx.max_wave
    minEm  = donEm.min_wave
    maxEm  = donEm.max_wave

    startingwave = int(max(minAcc, minEm))
    endingwave = int(min(maxAcc, maxEm))

    A = accEx.wave_value_pairs()
    D = donEm.wave_value_pairs()
    overlap = [(pow(wave, 4) * A[wave] * accEC * D[wave] / donCum) for
               wave in range(startingwave, endingwave+1)]

    return sum(overlap)


def forsterDist(donor, acceptor, n=1.4, k=2./3.):
    overlap = calculate_spectral_overlap(donor, acceptor)
    return .2108*(pow((k)*(pow(n, -4)*overlap), (1./6.)))


def fretEfficiency(distance, forster):
    return 1/(1+pow(distance/forster, 6))


def find_best_forster():
    from .models import Protein
    qs = Protein.objects.with_spectra()
    out = []
    withSpectra = []
    for p in qs:
        try:
            p.default_state.em_spectra.data
        except Exception:
            continue
        withSpectra.append(p)
    for donor in withSpectra:
        for acceptor in withSpectra:
            try:
                out.append([forsterDist(donor, acceptor), donor.name, acceptor.name])
            except Exception:
                continue
    return sorted(out, key=lambda x: x[0].real)

import math

import numpy as np


def halflife2decay(halflife):
    """convert exponential halflife to decay rate"""
    return math.log(2) / halflife


def decay_to_halflife(decay):
    """convert decay rate to exponential halflife"""
    return math.log(2) / decay


def timeconst_to_decay(timeconst):
    """convert exponential time constant to decay"""
    return 1 / timeconst


def decay_to_timeconst(decay):
    """convert exponential decay to time constant"""
    return 1 / decay


def observation_to_decay(observation, time, initial=1):
    """from observed fraction at timepoint, calculate exponential decay"""
    return math.log(observation / initial) / -time


def observation_to_halflife(observation, time, initial=1):
    """from observed fraction at timepoint, calculate exponential decay"""
    return math.log(2) / (math.log(observation / initial) / -time)


def ec_to_cross_section(EC):
    """convert from molar extinction coefficient to absorption cross section

    note: EC should be expressed in 1/(M*cm) and is the EC at the employed
    excitation wavelength. cross section returned is in units cm^2/molecule

    doi: 10.1038/ncomms1738

    Annu. Rev. Biophys. Biophys. Chern. 1989. 18:271-308

    Valeur, B., and M.N. Berberan-Santos. 2012. Molecular Fluorescence.
    Wiley-VCH Verlag GmbH & Co. KGaA, Weinheim, Germany.
    """
    avogadro = 6.0221409e23
    return 1000 * math.log(10) * EC / avogadro


def power_density_to_photon_flux(power, wavelength):
    """convert power density, in W/cm^2, i.e. J/(s cm^2) to photon flux
    wavelength should be expressed in microns

    in photons / (s * cm^2)
    """
    c = 2.9979245e14  # speed of light in microns / s
    h = 6.62607004e-34  # Planck's constant in J * s
    E = h * c / wavelength  # J per photon of wavelength lambda
    return power / E


def total_photon_prediction(k_bleach, EC, QY, power, wavelength):
    abs_cross_section = ec_to_cross_section(EC)  # cm^2/molecule
    photon_flux = power_density_to_photon_flux(power, wavelength)  # photons/(s*cm^2)
    k_excitation = abs_cross_section * photon_flux  # photons/molecule/s
    k_radiative_decay = QY * k_excitation  # photons/molecule/s
    N_emmission = k_radiative_decay / k_bleach  # photons/molecule
    return N_emmission


def fit_exp_linear(t, y, offset=0):
    y = y - offset
    y = np.log(y)
    K, A_log = np.polyfit(t, y, 1)
    A = np.exp(A_log)
    return A, K


def fit_points(t, y, offset=0):
    t = np.array(t)
    y = np.array(y)
    # Linear Fit (Note that we have to provide the y-offset ("C") value!!
    A, K = fit_exp_linear(t, y, offset)
    return A, -K

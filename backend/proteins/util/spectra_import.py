import numpy as np
from django.contrib.auth import get_user_model

from ..forms import SpectrumForm
from .helpers import zip_wave_data

User = get_user_model()


def get_stype(header):
    if "(2p)" in header.lower():
        return "2p"
    if any(x in header.lower() for x in ("(ab", "absorption")):
        return "ab"
    if any(x in header.lower() for x in ("(em", "emission")):
        return "em"
    if any(x in header.lower() for x in ("(ex", "excitation")):
        return "ex"


# FIXME: I kind of hate this function...
def import_spectral_data(waves, data, headers=None, categories=(), stypes=(), owner=None, minmax=None):
    """
    Take a vector of waves and a matrix of data, and import into database

    if cat or stype are strings, they will be assigned to all spectra
    they can also be a list of strings with the same length as the number of
    data columns (not including waves) in the file
    """
    if isinstance(categories, list | tuple) and len(categories):
        assert len(categories) == len(data), "provided category list not the same length as data"
    elif isinstance(categories, str):
        categories = [categories] * len(data)
    else:
        raise ValueError("Must provide category, or list of categories with same length as data")

    if isinstance(stypes, list | tuple) and len(stypes):
        assert len(stypes) == len(data), "provided subtypes list not the same length as data"
    elif isinstance(stypes, str):
        stypes = [stypes] * len(data)
    else:
        stypes = [None] * len(data)

    if not headers:
        headers = [None] * len(data)

    new_objects = []
    errors = []
    for datum, header, cat, stype in zip(data, headers, categories, stypes):
        if not (any(datum)) or all(np.isnan(datum)):
            print(f"skipping col {header} ... no data")
            continue

        if not stype:
            stype = get_stype(header)

        iowner = owner
        if not iowner:
            iowner = header.split(" (")[0]

        d = zip_wave_data(waves, datum, minmax)
        sf = SpectrumForm(
            {
                "data": d,
                "category": cat,
                "subtype": stype,
                "owner": iowner,
                "confirmation": True,
            }
        )
        if sf.is_valid():
            newob = sf.save(commit=False)
            newob.created_by = User.objects.first()
            newob.save()
            newob.owner.created_by = User.objects.first()
            newob.owner.save()
            new_objects.append(newob)
            print(f"Successfully imported {iowner}, {cat}, {stype}")
        else:
            errors.append((iowner, sf.errors))

    return new_objects, errors

import contextlib
import os
import re
from io import StringIO

import pandas as pd
import requests
from django.core.validators import URLValidator
from django.template.defaultfilters import slugify

from ..models import Filter
from .helpers import zip_wave_data

############################################
#       Importing Tools
############################################


BASEDIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


def reimport_filter_spectrum(obj):
    if obj.manufacturer.lower() == "semrock":
        text = fetch_semrock_part(obj.part)
    elif obj.manufacturer.lower() == "chroma":
        text = fetch_chroma_part(obj.part)
    waves, data, headers = text_to_spectra(text)
    D = zip_wave_data(waves, data[0])
    obj.spectrum.data = D
    obj.spectrum.save()


def add_filter_to_database(brand, part, user=None):
    if brand.lower() == "chroma":
        importer = import_chroma_spectra
    elif brand.lower() == "semrock":
        part = normalize_semrock_part(part)
        importer = import_semrock_spectra
    else:
        raise ValueError("unknown brand")

    with contextlib.suppress(Filter.DoesNotExist):
        f = Filter.objects.get(slug=slugify(brand + " " + part))
        return [f.spectrum], []

    new_objects, errors = importer(part=part)
    if new_objects and user:
        for spectrum in new_objects:
            spectrum.owner.created_by = user
            spectrum.owner.save()
    return new_objects, errors


def fetch_chroma_url(url):
    """
    url is either a Chroma ASCII text url:
    https://www.chroma.com/file/69074/download?token=nMjkwb45
    or a chroma token:
    'nMjkwb45'
    """

    urlv = URLValidator()

    try:
        urlv(url)
    except Exception as e:
        raise ValueError("invalid url for Chroma download") from e

    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise ValueError("ASCII download failed")


def check_chroma_for_part(part):
    part = part.replace("/", "-")
    chroma_url = "https://www.chroma.com/products/parts/"
    return requests.head(chroma_url + slugify(part))


def check_semrock_for_part(part):
    part = normalize_semrock_part(part)
    part = part.replace("/", "_").upper()
    semrock_url = "https://www.semrock.com/_ProductData/Spectra/"
    return requests.head(semrock_url + slugify(part) + "_Spectrum.txt")


def fetch_chroma_part(part):
    """Retrieve ASCII spectra for a chroma part number

    part is a string:
    'ET525/50m' or 'et525-50m'
    must resolve to a url such as:
    https://www.chroma.com/products/parts/et525-50m
    """
    from html.parser import HTMLParser

    class ChromaParser(HTMLParser):
        ready = 0
        url = False

        def handle_starttag(self, tag, attrs):
            if self.ready == 0:
                if tag == "span" and len(attrs):
                    for k, v in attrs:
                        if k == "class" and v == "file":
                            self.ready = 1
            elif self.ready == 1:
                if tag == "a" and len(attrs):
                    for k, v in attrs:
                        if k == "href" and "download?token" in v:
                            self.url = v
                            self.ready = 2

    part = part.replace("/", "-")
    chroma_url = "https://www.chroma.com/products/parts/"
    response = requests.get(chroma_url + slugify(part))
    if response.status_code != 200:
        raise ValueError(f"Could not retrieve Chroma part: {slugify(part)}")

    parser = ChromaParser()
    parser.feed(response.text)
    if parser.url:
        return fetch_chroma_url(parser.url)
    else:
        raise ValueError(f"Found Chroma part {slugify(part)}, but could not find file to download")


def normalize_semrock_part(part):
    return re.sub("-(25|35)$", "", part).upper()


def fetch_semrock_part(part):
    """Retrieve ASCII spectra for a semrock part number

     part is a string:
    'FF01-571/72' or  'FF01-571/72-25' (-25) will be clipped
     must resolve to a url such as:
     https://www.semrock.com/_ProductData/Spectra/FF01-571_72_Spectrum.txt
    """
    response = requests.get("https://www.semrock.com/FilterDetails.aspx?id=" + part)
    if response.status_code != 200:
        raise ValueError(f"Semrock part not valid for URL: {part}")

    try:
        url = (
            "https://www.semrock.com"
            + (str(response.content).split('" title="Click to Download ASCII')[0].split('href="')[-1])
        )
    except Exception as e:
        raise ValueError(f"Could not parse page for semrock part: {part}") from e

    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError(f"Could not retrieve data for Semrock part: {part}")

    txt = response.text
    if txt.startswith("Typical") and "Data format" in txt:
        txt = txt.split("Data format")[1]
        txt = "".join(txt.split("\n")[1:])
        txt = txt.split("---")[0]
        txt = txt.strip("\r").replace("\r", "\n")
        if len({len(x.split("\t")) for x in txt.split("\n")}) > 1:
            return "\n".join([x for x in txt.split("\n") if len(x.split("\t")) == 2])
    return txt


def all_numbers(df):
    return all(pd.api.types.is_numeric_dtype(t) for t in df.dtypes)


def read_csv_text(text) -> pd.DataFrame:
    fmt = {"sep": ";", "thousands": ".", "decimal": ","}
    df = pd.read_csv(StringIO(text), **fmt)
    if df.ndim != 2 or df.shape[1] < 2 or not all_numbers(df):
        fmt = {"sep": ",", "thousands": ",", "decimal": "."}
        df = pd.read_csv(StringIO(text), **fmt)
    if df.ndim != 2 or df.shape[1] < 2 or not all_numbers(df):
        raise ValueError("Could not parse text as valid spectra csv.")
    try:
        return pd.read_csv(StringIO(text), **fmt, header=None, dtype="float")
    except ValueError:
        return pd.read_csv(StringIO(text), **fmt, dtype="float")


def text_to_spectra(text, wavecol=0):
    """Convert text string into spectral data

    Args:
        text (str): string containing csv (usually) data
        wavecol (:obj:`int`, optional): column in csv that contains wavelengths

    Returns:
        tuple: (waves, outdata, headers).  waves is 1D, outdata is MxN, where M
            is the number of data columns and N is the number of wavelenghts.
            headers is 1D of length M, containing titles of data colums
    """
    df = read_csv_text(text)
    waves = df.iloc[:, wavecol].to_numpy(dtype="f")
    headers = [str(h) for i, h in enumerate(df.columns) if i != wavecol]
    outdata = df.drop(df.columns[wavecol], axis=1).to_numpy(dtype="f")
    return waves, outdata.T.tolist(), headers


def import_chroma_spectra(part=None, url=None, **kwargs):
    from ..util.spectra_import import import_spectral_data

    if isinstance(part, str):
        text = fetch_chroma_part(part)
        kwargs["owner"] = "Chroma " + part
        if part.lower().endswith(("m", "em")):
            kwargs["stypes"] = "bm"
        elif part.lower().endswith(("ex", "x", "bp")):
            kwargs["stypes"] = "bx"
        elif ("pc" in part.lower()) or ("bs" in part.lower()) or ("dc" in part.lower()):
            kwargs["stypes"] = "bs"
        elif "lp" in part.lower():
            kwargs["stypes"] = "lp"
        elif "sp" in part.lower():
            kwargs["stypes"] = "sp"
        if "stypes" not in kwargs:
            raise ValueError(f"Could not guess filter type for part {part}")
    elif url:
        if "owner" not in kwargs:
            raise ValueError('must provide argument "owner" when importing from url')
        if "stypes" not in kwargs:
            raise ValueError('must provide argument "stypes" when importing from url')
        text = fetch_chroma_url(url)
    else:
        ValueError("did not receive appropriate input to import_chroma_spectra")
    waves, data, headers = text_to_spectra(text)

    kwargs["categories"] = "f"
    newObjects, errors = import_spectral_data(waves, data, headers, **kwargs)
    for obj in newObjects:
        obj.source = "Chroma website"
        obj.save()
        obj.owner.manufacturer = "Chroma"
        obj.owner.part = part
        obj.owner.save()
    return newObjects, errors


def import_semrock_spectra(part=None, **kwargs):
    from ..util.spectra_import import import_spectral_data

    if isinstance(part, str):
        part = normalize_semrock_part(part)
        text = fetch_semrock_part(part)
        kwargs["owner"] = "Semrock " + part

        kwargs["stypes"] = "bp"  # default to bandpass
        if "sp" in part.lower():
            kwargs["stypes"] = "sp"
        if "lp" in part.lower():
            kwargs["stypes"] = "lp"
        if ("di" in part.lower()) or len([i for i in part if i == "/"]) > 1:
            kwargs["stypes"] = "bs"
        # if re.search(r'(\d+)/(\d+)', part):
        #     w1, w2 = re.search(r'(\d+)/(\d+)', part).groups()
        #     if int(w2) > int(w1):  # likely a dichroic
        #         kwargs['stypes'] = 'bs'
        if "stypes" not in kwargs:
            raise ValueError(f"Could not guess filter type for part {part}")
    else:
        ValueError("did not receive appropriate input to import_semrock_spectra")
    waves, data, headers = text_to_spectra(text)

    kwargs["categories"] = "f"

    new_objects, errors = import_spectral_data(waves, data, headers, **kwargs)
    for obj in new_objects:
        obj.source = "Semrock website"
        obj.save()
        obj.owner.manufacturer = "Semrock"
        obj.owner.part = part
        obj.owner.save()

    return new_objects, errors


#########################
#  Conversion
#########################

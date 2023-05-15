import contextlib
import io
import logging
import os
import re
from collections import Counter, OrderedDict
from math import isnan
from uuid import uuid4

import matplotlib.ticker as ticker
from django.core.cache import cache
from django.urls import reverse
from django.utils.safestring import mark_safe
from django.utils.text import slugify
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

logger = logging.getLogger(__name__)


@contextlib.contextmanager
def remember_cwd():
    curdir = os.getcwd()
    try:
        yield
    finally:
        os.chdir(curdir)


def create_slug_dict():
    from proteins.models import Protein

    slugs = OrderedDict(Protein.objects.all().values_list("name", "slug"))
    for item in Protein.objects.exclude(aliases=[]).values_list("aliases", "slug"):
        if item[0]:
            for alias in item[0]:
                slugs.update({alias: item[1]})
    return OrderedDict(sorted(slugs.items(), key=lambda x: len(x[0]), reverse=True))


def link_excerpts(excerpts_qs, obj_name=None, aliases=()):
    if not excerpts_qs:
        return None
    excerpt_list = list(excerpts_qs)
    slug_dict = cache.get_or_set("slug_dict", create_slug_dict, 60)
    for excerpt in excerpt_list:
        for name in slug_dict or {}:
            if len(name) <= 1:
                continue
            if name == obj_name or (aliases and name in aliases):
                excerpt.content = mark_safe(
                    re.sub(
                        rf"(?<=[\s(])(?<!>){name}(?!.\d)(?!<)",
                        f"<strong>{name}</strong>",
                        excerpt.content,
                    )
                )
            else:
                excerpt.content = mark_safe(
                    re.sub(
                        rf"(?<=[\s(])(?<!>){name}(?!.\d)(?!<)",
                        '<a href="{}" class="text-info">{}</a>'.format(
                            reverse("proteins:protein-detail", args=[slug_dict[name]]),
                            name,
                        ),
                        excerpt.content,
                    )
                )
    return excerpt_list


def most_favorited(max_results=20):
    from favit.models import Favorite
    from proteins.models import Protein

    qs = Favorite.objects.for_model(Protein)
    fave_counts = Counter(qs.values_list("target_object_id", flat=True))
    fave_items = dict(fave_counts.most_common(max_results))
    qs = Protein.objects.filter(id__in=fave_items.keys()).values("id", "name", "slug")
    D = {q.pop("id"): q for q in qs}

    od = OrderedDict()
    for prot_id, count in fave_items.items():
        od[prot_id] = D[prot_id]
        od[prot_id]["count"] = count
    return od


def merge_proteins(merge_prot, into_prot):
    # need to update favorites
    from favit.models import Favorite

    for fav in Favorite.objects.filter(target_object_id=merge_prot.id):
        try:
            fav.target_object_id = into_prot.id
            fav.save()
        except Exception:
            fav.delete()
    into_prot.aliases.append(merge_prot.name)
    # need a LOT more work...


def getprot(protein_name, visible=False):
    from proteins.models import Protein

    # assume that the slug is always the slugified name
    try:
        if visible:
            return Protein.visible.get(slug=slugify(protein_name))
        else:
            return Protein.objects.get(slug=slugify(protein_name))
    except Exception:
        pass
    return Protein.objects.get(aliases__contains=[protein_name])


def getmut(protname2, protname1=None, ref=None):
    from fpseq.mutations import get_mutations

    a = getprot(protname2)
    if protname1:
        b = getprot(protname1)
    else:
        b = a.lineage.parent.protein
    ref = getprot(ref).seq if ref else None
    return get_mutations(b.seq, a.seq, ref)


def showalign(protname2, protname1=None):
    a = getprot(protname2)
    if protname1:
        b = getprot(protname1)
    else:
        b = a.lineage.parent.protein
    print(b.seq.align_to(a.seq))


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
    minmax = minmax or (150, 1800)
    return [list(i) for i in zip(waves, data) if (minmax[0] <= i[0] <= minmax[1]) and not isnan(i[1])]


def wave_to_hex(wavelength, gamma=1):
    """This converts a given wavelength into an approximate RGB value.
    The given wavelength is in nanometers.
    The range of wavelength is 380 nm through 750 nm.

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    """
    if not wavelength:
        return "#000"

    wavelength = float(wavelength)
    if 520 <= wavelength:
        wavelength += 40

    if wavelength < 380:
        r = 0.05
        g = 0.0
        b = 0.15
    elif wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        r = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        g = 0.0
        b = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        r = 0.0
        g = ((wavelength - 440) / (490 - 440)) ** gamma
        b = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        r = 0.0
        g = 1.0
        b = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        r = ((wavelength - 510) / (580 - 510)) ** gamma
        g = 1.0
        b = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        r = 1.0
        g = (-(wavelength - 645) / (645 - 580)) ** gamma
        b = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (770 - wavelength) / (770 - 645)
        r = (1.0 * attenuation) ** gamma
        g = 0.0
        b = 0.0
    else:
        r = 0.18
        g = 0.0
        b = 0.05
    r *= 255
    g *= 255
    b *= 255
    return f"#{int(r):02x}{int(g):02x}{int(b):02x}"


# def wave_to_hex(wave):
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
    if re.search("^m[A-Z]", name):
        return name.lstrip("m")
    if name.lower().startswith("monomeric"):
        name = name[9:]
    return name.lstrip(" ")


def get_base_name(name):
    '''return core name of protein, stripping prefixes like "m" or "Tag"'''

    # remove PA/(Pa), PS, PC, from beginning
    if name.startswith(("PA", "Pa", "PS", "Ps", "PC", "pc", "rs")):
        name = name[2:]

    if re.match("LSS", name):
        name = name[3:].lstrip("-")

    # remove m (if next letter is caps) or monomeric
    if re.match("m[A-Z]", name):
        name = name[1:]

    # get rid of Td or td
    if re.match("[Tt][Dd][A-Z]", name):
        name = name[2:]

    if name.lower().startswith("monomeric"):
        name = name[9:]
    if name.lower().startswith("tag"):
        name = name[3:]

    # remove E at beginning (if second letter is caps)
    if re.match("E[A-Z]", name):
        name = name[1:]
    # remove S at beginning (if second letter is caps)
    if re.match("S[A-Z]", name):
        name = name[1:]
    # remove T- at beginning (if second letter is caps)
    if re.match("T-", name):
        name = name[2:]

    name = name.lstrip("-").lstrip(" ")

    return name


# ###########################################
#       Spectral Functions
# ###########################################


def calculate_spectral_overlap(donor, acceptor):
    accEx = acceptor.default_state.ex_spectrum
    accEC = acceptor.default_state.ext_coeff
    donEm = donor.default_state.em_spectrum
    # donQY  = donor.default_state.qy
    donCum = sum(donEm.y)

    minAcc = accEx.min_wave
    maxAcc = accEx.max_wave
    minEm = donEm.min_wave
    maxEm = donEm.max_wave

    startingwave = int(max(minAcc, minEm))
    endingwave = int(min(maxAcc, maxEm))

    A = accEx.wave_value_pairs()
    D = donEm.wave_value_pairs()
    overlap = [(pow(wave, 4) * A[wave] * accEC * D[wave] / donCum) for wave in range(startingwave, endingwave + 1)]

    return sum(overlap)


def forsterDist(donor, acceptor, n=1.329, k=2 / 3):
    overlap = calculate_spectral_overlap(donor, acceptor)
    return (
        overlap * 1e-15,
        0.2108 * pow(donor.default_state.qy * k * pow(n, -4) * overlap, (1 / 6)),
    )


def fretEfficiency(distance, forster):
    return 1 / (1 + pow(distance / forster, 6))


def forster_list():
    from ..models import Protein

    qs = (
        Protein.objects.with_spectra()
        .filter(agg=Protein.MONOMER, switch_type=Protein.BASIC)
        .select_related("default_state")
        .prefetch_related("default_state__spectra")
    )
    out = []
    withSpectra = []
    for p in qs:
        try:
            _ = p.default_state.em_spectrum.data
        except Exception:
            continue
        withSpectra.append(p)
    for donor in withSpectra:
        for acceptor in withSpectra:
            try:
                if (
                    (acceptor.default_state.ex_max > donor.default_state.ex_max)
                    and acceptor.default_state.ext_coeff
                    and donor.default_state.qy
                ):
                    overlap, r0 = forsterDist(donor, acceptor)
                    out.append(
                        {
                            "donor": "<a href='{}'>{}{}</a>".format(
                                donor.get_absolute_url(),
                                donor.name,
                                f"<sub>{donor.cofactor.upper()}</sub>" if donor.cofactor else "",
                            ),
                            "acceptor": "<a href='{}'>{}{}</a>".format(
                                acceptor.get_absolute_url(),
                                acceptor.name,
                                f"<sub>{acceptor.cofactor.upper()}</sub>" if acceptor.cofactor else "",
                            ),
                            "donorPeak": donor.default_state.ex_max,
                            "acceptorPeak": acceptor.default_state.ex_max,
                            "emdist": acceptor.default_state.em_max - donor.default_state.em_max,
                            "donorQY": donor.default_state.qy,
                            "acceptorQY": acceptor.default_state.qy,
                            "acceptorEC": f"{acceptor.default_state.ext_coeff:,}",
                            "overlap": round(overlap, 2),
                            "forster": round(r0.real, 2),
                            "forsterQYA": round(r0.real * acceptor.default_state.qy, 2),
                        }
                    )
            except Exception:
                continue
    return sorted(out, key=lambda x: x["forster"], reverse=True)


def spectra_fig(
    spectra,
    format="svg",
    output=None,
    xlabels=True,
    ylabels=False,
    xlim=None,
    fill=True,
    transparent=True,
    grid=False,
    title=False,
    info=None,
    figsize=(12, 3),
    **kwargs,
):
    if not spectra:
        return None

    alph = kwargs.pop("alpha", None)
    colr = kwargs.pop("color", None)
    twitter = kwargs.pop("twitter", 0)

    logger.debug("spectra_fig called on {}".format(",".join([str(spec.id) for spec in spectra])))
    if twitter:
        xlabels = False
        transparent = False
        figsize = (12, 6)
        xlim = (400, 760)
        xlabels = True
    fig = Figure(figsize=figsize, dpi=70)
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    if transparent:
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)
    if not xlim:
        xlim = (min([s.min_wave for s in spectra]), max([s.max_wave for s in spectra]))
    for spec in spectra:
        color = spec.color() if not colr else colr
        if fill:
            alpha = 0.5 if not alph else float(alph)
            ax.fill_between(*list(zip(*spec.data)), color=color, alpha=alpha, url="http://google.com=", **kwargs)
        else:
            alpha = 1 if not alph else float(alph)
            ax.plot(*list(zip(*spec.data)), alpha=alpha, color=spec.color(), **kwargs)
    ax.set_ylim((-0.005, 1.025))
    ax.set_xlim(xlim)
    if twitter:
        ax.set_ylim((0, 1.07))
    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    if grid:
        ax.grid(color="gray", axis="both", alpha=0.15, which="both")
        ax.xaxis.grid(True, which="minor")
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        ax.set_axisbelow(True)

    pos = [0, 0.017, 0.97, 0.98]
    if xlabels:
        ax.spines["bottom"].set_linewidth(0.4)
        ax.spines["bottom"].set_color((0.5, 0.5, 0.5))
        ax.tick_params(axis="x", colors=(0.2, 0.2, 0.2), length=0)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
        pos[0] = 0.02
        pos[1] = 0.08
        pos[3] -= 0.065
    else:
        ax.spines["bottom"].set_visible(False)
        ax.get_xaxis().set_ticks([])
    if ylabels:
        ax.tick_params(axis="y", colors=(0.2, 0.2, 0.2), length=0)
        ax.spines["left"].set_linewidth(0.4)
        ax.spines["left"].set_color((0.5, 0.5, 0.5))
        pos[0] = 0.025
        pos[2] = 0.96
        pos[3] -= 0.01
    else:
        ax.spines["left"].set_visible(False)
        ax.get_yaxis().set_ticks([])
        pos[0] = 0.015

    ax.set_position(pos)
    if title:
        font = {
            "family": "sans-serif",
            "color": "black",
            "weight": "normal",
            "size": 18,
        }
        ax.text(xlim[0] + 2, 0.97, title, va="top", fontdict=font, alpha=0.5)
        if info:
            font["size"] = 14
            ax.text(xlim[0] + 2, 0.85, info, va="top", fontdict=font, alpha=0.5)

    if not output:
        output = io.BytesIO()
    canvas.print_figure(output, format=format)
    return output


def wipe_bad_uuids():
    """get rid of old uuids in version histories"""
    import json

    from reversion.models import Version

    for version in Version.objects.all():
        data = json.loads(version.serialized_data)
        for item in data:
            if "fields" in item and "uuid" in item["fields"] and len(item["fields"]["uuid"]) > 5:
                item["fields"]["uuid"] = "-----"
        version.serialized_data = json.dumps(data)
        version.save()

import random
from typing import TYPE_CHECKING, cast

import factory
import factory.fuzzy
import numpy as np
from django.utils.text import slugify

from fpseq import FPSeq
from proteins.util.helpers import wave_to_hex

from .models import Filter, FilterPlacement, Microscope, OpticalConfig, Protein, Spectrum, State

if TYPE_CHECKING:
    import factory.builder

# fmt: off
_BP_TYPE = {Filter.BP, Filter.BPM, Filter.BPX}

FRUITS = [
    'Apple', 'Banana', 'Orange', 'Strawberry', 'Mango', 'Grapes', 'Pineapple', 'Watermelon', 'Kiwi',
    'Pear', 'Cherry', 'Peach', 'Plum', 'Lemon', 'Lime', 'Blueberry', 'Raspberry', 'Blackberry',
    'Pomegranate', 'Coconut', 'Avocado', 'Grapefruit', 'Cantaloupe', 'Fig', 'Guava', 'Honeydew',
    'Lychee', 'Mandarin', 'Nectarine', 'Passionfruit', 'Papaya', 'Apricot', 'Cranberry', 'Dragonfruit',
    'Starfruit', 'Persimmon', 'Rambutan', 'Tangerine', 'Clementine', 'Mulberry', 'Cactus',
    'Quince', 'Date', 'Jackfruit', 'Kumquat', 'Lingonberry', 'Loquat', 'Mangosteen', 'Pitaya',
    'Plantain', 'PricklyPear', 'Tamarind', 'UgliFruit', 'Yuzu', 'Boysenberry', 'Cherimoya',
]
LIQUERS = [
    'Campari', 'Chartreuse', 'Curacao', 'Midori', 'Cointreau', 'Galliano',
    'Sambuca', 'Cassis', 'Pernod', 'Amaro', 'TripleSec', 'Frangelico',
    'Baileys', 'Kahlua', 'Chambord', 'Cynar', 'Aperol', 'StGermain',
    'Benedictine', 'Drambuie',
]
AMINO_ACIDS = 'ARNDCEQGHILKMFPSTWYV'
COMMON_ORGANISMS = [
    (6100, 'Aequorea victoria', 'hydrozoans'),
    (86600, 'Discosoma sp.', 'coral anemones'),
    (6118, 'Entacmaea quadricolor', 'sea anemones'),
    (46758, 'Lobophyllia hemprichii', 'stony corals'),
    (86521, 'Clavularia sp.', 'soft corals'),
    (32630, 'synthetic construct', 'other sequences'),
    (63558, 'Montastraea cavernosa', 'stony corals'),
    (301887, 'Echinophyllia sp. SC22', 'stony corals'),
    (1076, 'Rhodopseudomonas palustris', 'a-proteobacteria'),
    (147615, 'Aequorea macrodactyla', 'hydrozoans'),
    (45264, 'Acropora millepora', 'stony corals'),
    (496660, 'Verrillofungia concinna', 'stony corals'),
    (3702, 'Arabidopsis thaliana', 'eudicots'),
    (105402, 'Zoanthus sp.', 'mat anemones'),
    (406427, 'Anthoathecata', 'hydrozoans'),
    (258594, 'Rhodopseudomonas palustris CGA009', 'a-proteobacteria'),
    (175771, 'Heteractis crispa', 'sea anemones'),
    (1495449, 'Olindias formosus', 'hydrozoans'),
    (321802, 'Montipora sp. 20', 'stony corals')
]
REAL_PDBS = [
    '1RRX', '2Z6X', '4B30', '3LVA', '6AA7', '2A46', '2A48', '2A47', '2C9I', '1XQM', '1EMA',
    '1RM9', '1RMP', '1RMM', '1H6R', '4DKM', '2FWQ', '3LF4', '6RHF', '3EVP', '2Q57', '2WSO',
    '7E70', '7E6Z', '7E6Y', '7E6X', '1HUY', '3DPW', '3SV5', '5WJ2', '3RWA', '2DD7', '1JBZ',
    '1JBY', '2VZX', '1ztu', '3ST3', '3ST2', '3ST4', '2POX', '2GX0', '2IE2', '2IOV', '2GX2',
]
# fmt: on

############### HELPERS ###############


def _protein_seq(length: int = 230) -> FPSeq:
    return FPSeq("M" + "".join(random.choices(AMINO_ACIDS, k=length - 1)))


def _skewed_gaussian(x, mean, sigma: float = 50, skewness: float = 0, amplitude: float = 1):
    t = (x - mean) / sigma
    y = 2 * np.exp(-0.5 * t**2) * (1 + np.tanh(skewness * t))
    y /= np.max(y)
    return amplitude * y


def _mock_spectrum(mean, sigma: float = 40, min_wave=300, max_wave=900, type="ex"):
    """Return a mock spectrum in the form [[x1, y1], [x2, y2], ...]."""
    match type:
        case "ex":
            skewness = -3.3
            mean += 17
        case "em":
            skewness = 3
            mean -= 21
            sigma = 50
        case _:
            skewness = 0
    x = np.arange(min_wave, max_wave + 1)
    y = _skewed_gaussian(x, mean, sigma, skewness)
    return np.stack([x, y], axis=1).tolist()


def _mock_bandpass(bandcenter, bandwidth, min_wave=300, max_wave=900, transmission=0.9):
    half_width = bandwidth / 2
    max_wave = max(max_wave, bandcenter + half_width)
    min_wave = min(min_wave, bandcenter - half_width)
    x = np.arange(min_wave, max_wave + 1)
    y = np.zeros_like(x, dtype="float")
    y[(bandcenter - half_width <= x) & (x <= bandcenter + half_width)] = transmission
    return np.stack([x, y], axis=1).tolist()


def _mock_edge_filter(edge, subtype, min_wave=300, max_wave=900, transmission=0.9):
    max_wave = max(max_wave, edge + 10)
    min_wave = min(min_wave, edge - 10)
    x = np.arange(min_wave, max_wave + 1)
    y = np.zeros_like(x, dtype="float")
    if subtype == Filter.SP:
        y[x < edge] = transmission
    else:
        y[x > edge] = transmission
    return np.stack([x, y], axis=1).tolist()


def _build_spectral_data(resolver: factory.builder.Resolver):
    subtype = getattr(resolver, "subtype", None)

    if (owner_state := getattr(resolver, "owner_state", None)) is not None:
        owner_state = cast(State, owner_state)
        if subtype == "ex":
            return _mock_spectrum(owner_state.ex_max, type="ex")
        elif subtype == "em":
            return _mock_spectrum(owner_state.em_max, type="em")

    if (owner_filter := getattr(resolver, "owner_filter", None)) is not None:
        owner_filter = cast(Filter, owner_filter)
        if subtype in _BP_TYPE:
            center = owner_filter.bandcenter or resolver.factory_parent.bandcenter
            width = owner_filter.bandwidth or resolver.factory_parent.bandwidth
            return _mock_bandpass(center, width)
        edge = owner_filter.edge or resolver.factory_parent.edge
        return _mock_edge_filter(edge, subtype)
    return _mock_spectrum(random.randint(450, 650), type="")


############### FACTORIES ###############


class OrganismFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = "proteins.Organism"
        django_get_or_create = ("id",)

    id = factory.Iterator([o[0] for o in COMMON_ORGANISMS])
    scientific_name = factory.Iterator([o[1] for o in COMMON_ORGANISMS])
    division = factory.Iterator([o[2] for o in COMMON_ORGANISMS])


class SpectrumOwnerFactory(factory.django.DjangoModelFactory):
    class Meta:
        abstract = True
        django_get_or_create = ("name", "slug")

    name = factory.Sequence(lambda n: f"TestSpectrum{n}")
    slug = factory.LazyAttribute(lambda o: slugify(o.name.replace("/", "-")))


class FluorophoreFactory(SpectrumOwnerFactory):
    class Meta:
        abstract = True
        django_get_or_create = ("name", "slug")

    name = factory.Sequence(lambda n: f"TestFluorophore{n}")
    ex_max = factory.Faker("pyint", min_value=400, max_value=650)
    em_max = factory.LazyAttribute(lambda o: o.ex_max + random.randint(10, 100))
    ext_coeff = factory.Faker("pyint", min_value=10, max_value=300000)
    qy = factory.Faker("pyfloat", min_value=0.1, max_value=1.0)
    pka = factory.Faker("pyfloat", min_value=2.0, max_value=12.0)
    lifetime = factory.Faker("pyfloat", min_value=1, max_value=16)
    brightness = factory.LazyAttribute(lambda o: float(round(o.ext_coeff * o.qy / 1000, 2)))
    emhex = factory.LazyAttribute(lambda o: "#000" if o.is_dark else wave_to_hex(o.em_max))
    exhex = factory.LazyAttribute(lambda o: wave_to_hex(o.ex_max))
    is_dark = False


class StateFactory(FluorophoreFactory):
    class Meta:
        model = State
        django_get_or_create = ("name", "slug")

    name = "default"
    slug = factory.LazyAttribute(lambda o: f"{o.protein.slug}_{slugify(o.name)}")
    maturation = factory.Faker("pyfloat", min_value=0, max_value=1600)
    protein = factory.SubFactory("proteins.factories.ProteinFactory", default_state=None)

    ex_spectrum = factory.RelatedFactory(
        "proteins.factories.SpectrumFactory",
        factory_related_name="owner_state",
        subtype="ex",
        category="p",
    )
    em_spectrum = factory.RelatedFactory(
        "proteins.factories.SpectrumFactory",
        factory_related_name="owner_state",
        subtype="em",
        category="p",
    )


class DyeFactory(FluorophoreFactory):
    class Meta:
        model = "proteins.Dye"


class ProteinFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = Protein
        django_get_or_create = ("name", "slug")

    # name = factory.LazyAttribute(lambda o: f"{o.agg}{random.choice(FRUITS + LIQUERS)}")
    name = factory.Sequence(lambda n: f"mTest{n}")
    slug = factory.LazyAttribute(lambda o: slugify(o.name))
    seq = factory.LazyFunction(_protein_seq)
    seq_validated = factory.Faker("boolean", chance_of_getting_true=75)
    agg = factory.fuzzy.FuzzyChoice(Protein.AGG_CHOICES, getter=lambda c: c[0])
    pdb = factory.LazyFunction(lambda: random.choices(REAL_PDBS, k=random.randint(0, 2)))
    parent_organism = factory.SubFactory(OrganismFactory)
    primary_reference = factory.SubFactory("references.factories.ReferenceFactory")
    default_state = factory.RelatedFactory(StateFactory, factory_related_name="protein")


class SpectrumFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = Spectrum

    category = factory.fuzzy.FuzzyChoice(Spectrum.CATEGORIES, getter=lambda c: c[0])
    subtype = factory.LazyAttribute(lambda o: random.choice(Spectrum.category_subtypes[o.category]))
    data = factory.LazyAttribute(_build_spectral_data)


class MicroscopeFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = Microscope

    name = factory.Sequence(lambda n: f"TestMicroscope{n}")


class FilterFactory(SpectrumOwnerFactory):
    class Meta:
        model = Filter
        exclude = ("subtype",)

    name = factory.Sequence(lambda n: f"TestFilter{n}")
    subtype = factory.fuzzy.FuzzyChoice(Spectrum.category_subtypes["f"])
    bandcenter = factory.LazyAttribute(lambda o: random.randint(400, 900) if o.subtype in _BP_TYPE else None)
    bandwidth = factory.LazyAttribute(lambda o: random.randint(10, 20) if o.subtype in _BP_TYPE else None)
    edge = factory.LazyAttribute(lambda o: random.randint(400, 900) if o.subtype not in _BP_TYPE else None)

    spectrum = factory.RelatedFactory(
        "proteins.factories.SpectrumFactory",
        factory_related_name="owner_filter",
        category="f",
        subtype=factory.SelfAttribute("..subtype"),
    )


class OpticalConfigFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = OpticalConfig

    name = factory.Sequence(lambda n: f"TestOC{n}")
    microscope = factory.SubFactory(MicroscopeFactory)
    # light
    # camera
    # laser


# this is a through table
class FilterPlacementFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = FilterPlacement

    filter = factory.SubFactory(FilterFactory)
    config = factory.SubFactory(OpticalConfigFactory)
    path = factory.fuzzy.FuzzyChoice(FilterPlacement.PATH_CHOICES, getter=lambda c: c[0])


class OpticalConfigWithFiltersFactory(OpticalConfigFactory):
    ex_filter = factory.RelatedFactory(
        FilterPlacementFactory,
        factory_related_name="config",
        filter__subtype=Filter.BP,
        path=FilterPlacement.EX,
    )
    bs_filter = factory.RelatedFactory(
        FilterPlacementFactory,
        factory_related_name="config",
        filter__subtype=Filter.BS,
        path=FilterPlacement.BS,
    )
    em_filter = factory.RelatedFactory(
        FilterPlacementFactory,
        factory_related_name="config",
        filter__subtype=Filter.BP,
        path=FilterPlacement.EM,
    )

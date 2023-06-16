import random

import factory
import faker.providers

from .models.protein import Protein

# fmt: off
FRUITS: set[str] = {
    'Apple', 'Banana', 'Orange', 'Strawberry', 'Mango', 'Grapes', 'Pineapple', 'Watermelon', 'Kiwi',
    'Pear', 'Cherry', 'Peach', 'Plum', 'Lemon', 'Lime', 'Blueberry', 'Raspberry', 'Blackberry',
    'Pomegranate', 'Coconut', 'Avocado', 'Grapefruit', 'Cantaloupe', 'Fig', 'Guava', 'Honeydew',
    'Lychee', 'Mandarin', 'Nectarine', 'Passionfruit', 'Papaya', 'Apricot', 'Cranberry', 'Dragonfruit',
    'Starfruit', 'Persimmon', 'Rambutan', 'Tangerine', 'Clementine', 'Mulberry', 'Cactus',
    'Quince', 'Date', 'Jackfruit', 'Kumquat', 'Lingonberry', 'Loquat', 'Mangosteen', 'Pitaya',
    'Plantain', 'PricklyPear', 'Tamarind', 'UgliFruit', 'Yuzu', 'Boysenberry', 'Cherimoya',
}
LIQUERS = {
    'Campari', 'Chartreuse', 'Curacao', 'Midori', 'Cointreau', 'Galliano',
    'Sambuca', 'Cassis', 'Pernod', 'Amaro', 'TripleSec', 'Frangelico',
    'Baileys', 'Kahlua', 'Chambord', 'Cynar', 'Aperol', 'StGermain',
    'Benedictine', 'Drambuie',
}
AMINO_ACIDS = 'ARNDCEQGHILKMFPSTWYV'
# fmt: on

ORGANISMS = {
    287157: (287157, "Acropora aculeus", "stony corals"),
    70779: (70779, "Acropora digitifera", "stony corals"),
    526283: (526283, "Acropora eurystoma", "stony corals"),
    55974: (55974, "Acropora hyacinthus", "stony corals"),
    45264: (45264, "Acropora millepora", "stony corals"),
    70781: (70781, "Acropora nobilis", "stony corals"),
    140239: (140239, "Acropora pulchra", "stony corals"),
    258444: (258444, "Acropora sp. #30", "stony corals"),
    70783: (70783, "Acropora tenuis", "stony corals"),
    6106: (6106, "Actinia equina", "sea anemones"),
    1246302: (1246302, "Aequorea australis", "hydrozoans"),
    210840: (210840, "Aequorea coerulescens", "hydrozoans"),
    147615: (147615, "Aequorea macrodactyla", "hydrozoans"),
    6100: (6100, "Aequorea victoria", "hydrozoans"),
    89882: (89882, "Agaricia agaricites", "stony corals"),
    165097: (165097, "Agaricia fragilis", "stony corals"),
    105399: (105399, "Anemonia majano", "sea anemones"),
    6108: (6108, "Anemonia sulcata", "sea anemones"),
    7937: (7937, "Anguilla japonica", "bony fishes"),
    406427: (406427, "Anthoathecata", "hydrozoans"),
}


@factory.Faker.add_provider
class FPbaseProvider(faker.providers.BaseProvider):
    def fp_name(self) -> str:
        prefix: str = self.random_element({"m", "d", "td", ""})
        return prefix + self.random_element(FRUITS | LIQUERS)

    def fp_seq(self, length: int = 230) -> str:
        return "M" + "".join(self.random_elements(AMINO_ACIDS, length=length))

    def organism_id(self) -> int:
        return random.choice(list(ORGANISMS))


class OrganismFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = "proteins.Organism"

    id = factory.Faker("organism_id")
    scientific_name = factory.LazyAttribute(lambda o: ORGANISMS[o.id][1])
    division = factory.LazyAttribute(lambda o: ORGANISMS[o.id][2])


class ProteinFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = "proteins.Protein"

    name = factory.Faker("fp_name")
    seq = factory.Faker("fp_seq", length=230)
    agg = factory.LazyAttribute(lambda o: next((x for x in ("m", "d", "td") if o.name.startswith(x)), ""))
    parent_organism = factory.SubFactory(OrganismFactory)

    @classmethod
    def get_or_create(cls, **kwargs):
        try:
            return Protein.objects.get(**kwargs), False
        except cls._meta.model.DoesNotExist:
            return cls.create(**kwargs), True

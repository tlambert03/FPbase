import factory
import factory.fuzzy

REAL_DOIS = [
    "10.1016/0378-1119(92)90691-h",
    "10.1126/science.8303295",
    "10.1073/pnas.91.26.12501",
    "10.1038/nbt0295-151",
    "10.1038/373663b0",
    "10.1016/0378-1119(95)00685-0",
    "10.1016/0378-1119(95)00768-7",
    "10.1016/s0960-9822(02)00450-5",
    "10.1038/nbt0396-315",
    "10.1126/science.273.5280.1392",
    "10.1038/nbt1096-1246",
    "10.1006/bbrc.1996.1573",
    "10.1093/nar/24.22.4592",
    "10.1099/00221287-143-2-303",
    "10.1038/nsb0597-361",
    "10.1021/bi970563w",
    "10.1016/s0006-3495(97)78307-3",
    "10.1074/jbc.272.45.28545",
    "10.1016/s0091-679x(08)61946-9",
    "10.1074/jbc.273.14.8212",
    "10.1146/annurev.biochem.67.1.509",
    "10.1038/28190",
    "10.1016/s0969-2126(98)00127-0",
    "10.2144/98255bt01",
    "10.1073/pnas.96.5.2135",
    "10.1109/5.771073",
    "10.1038/10904",
    "10.1038/13657",
    "10.1016/0014-5793(95)00557-p",
    "10.1093/nar/28.16.e78",
    "10.1073/pnas.97.22.11984",
    "10.1073/pnas.97.22.11990",
    "10.1126/science.290.5496.1585",
    "10.1038/81992",
    "10.1073/pnas.97.26.14091",
    "10.1073/pnas.051636098",
    "10.1074/jbc.m102815200",
    "10.1016/s0014-5793(01)02930-1",
    "10.1093/emboj/20.21.5853",
    "10.1016/s0301-0104(01)00486-4",
    "10.1016/0014-5793(94)80472-9",
    "10.1038/nbt0102-83",
    "10.1038/nbt0102-87",
    "10.1007/s10126-001-0081-7",
    "10.1016/s1096-4959(02)00025-8",
    "10.1073/pnas.062552299",
    "10.1126/science.1068539",
    "10.1073/pnas.082243699",
    "10.1073/pnas.182157199",
    "10.1126/science.1074952",
    "10.1073/pnas.202320599",
    "10.1074/jbc.m209524200",
    "10.1042/bj20021191",
    "10.1021/bi026609p",
    "10.1016/s0896-6273(02)01099-1",
    "10.1186/1472-6750-3-5",
]


# This is not ideal, the Reference.save() method will hit the actual REST API
# but after trying a bit to override the save() method specifically for this
# factory, I gave up and just used the real DOI's.
class ReferenceFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = "references.Reference"
        django_get_or_create = ("doi",)

    doi = factory.fuzzy.FuzzyChoice(REAL_DOIS)


class AuthorFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = "references.Author"
        django_get_or_create = ("family", "initials")

    family = factory.Faker("last_name")
    given = factory.Faker("first_name")
    initials = factory.LazyAttribute(lambda o: o.given[0])

# from django.utils.translation import gettext as _
import ast
import re

from Bio import Data, Seq
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator

from fpseq.mutations import Mutation

from .fields import Spectrum

validate_doi = RegexValidator(r"^10.\d{4,9}/[-._;()/:a-zA-Z0-9]+$", "Not a valid DOI string")
validate_uniprot = RegexValidator(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",
    "Not a valid UniProt Accession",
)

UNAMBIGUOUS_DNA_LETTERS = Seq.Seq("GATC")
IUPAC_PROTEIN_LETTERS = Seq.Seq("ACDEFGHIKLMNPQRSTVWY")


def validate_mutationset(mutset):
    errors = []
    for mut in re.split(" |,|/|\\\\", mutset):
        if mut:
            try:
                validate_mutation(mut)
            except ValidationError as e:
                errors.append(ValidationError(f"Bad Mutation String: {e}", code="badmut"))
    if errors:
        raise ValidationError(errors)


def validate_mutation(code):
    try:
        code = code.strip()
        m = Mutation.from_str(code)
        if str(m) != code:
            raise ValueError(f"Parsed mutation ({m}) different than input ({code})")
    except ValueError as e:
        raise ValidationError(f"Invalid mutation: {e}") from e


def cdna_sequence_validator(seq):
    badletters = [letter for letter in seq if letter not in UNAMBIGUOUS_DNA_LETTERS]
    if len(badletters):
        raise ValidationError(f"Invalid DNA letters: {''.join(set(badletters))}")
    try:
        Seq.Seq(seq, UNAMBIGUOUS_DNA_LETTERS).translate()
    except Data.CodonTable.TranslationError as e:
        raise ValidationError(e) from e


def protein_sequence_validator(seq):
    seq = "".join(str(seq).split()).upper()  # remove whitespace
    badletters = [letter for letter in seq if letter not in IUPAC_PROTEIN_LETTERS]
    if len(badletters):
        badletters = set(badletters)
        raise ValidationError(f"Invalid letter(s) found in amino acid sequence: {''.join(badletters)}")


def validate_spectrum(value):
    if not value:
        return None
    if isinstance(value, Spectrum):
        return
    try:
        obj = ast.literal_eval(value)
    except Exception as e:
        raise ValidationError("Invalid input for a Spectrum instance") from e
    if not isinstance(obj, list):  # must be a list
        raise ValidationError("Spectrum object must be of type List")
    if not all(isinstance(elem, list | tuple) for elem in obj):  # must be list of lists
        raise ValidationError("Spectrum object must be a list of lists or tuples")
    for elem in obj:
        if len(elem) != 2:
            raise ValidationError("All elements in Spectrum list must have two items")
        if not all(isinstance(n, int | float) for n in elem):
            raise ValidationError("All items in Spectrum list elements must be numbers")

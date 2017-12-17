from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from proteins.models import Spectrum
import ast

validate_doi = RegexValidator(r"^10.\d{4,9}/[-._;()/:a-zA-Z0-9]+$", 'Not a valid DOI string')


def validate_spectrum(value):
    if not value:
        return None
    if isinstance(value, Spectrum):
        return
    try:
        obj = ast.literal_eval(value)
    except Exception:
        raise ValidationError("Invalid input for a Spectrum instance")
    if not isinstance(obj, list):                           # must be a list
        raise ValidationError("Spectrum object must be of type List")
    if not all(isinstance(elem, list) for elem in obj):      # must be list of lists
        raise ValidationError("Spectrum object must be a list of lists")
    for elem in obj:
        if not len(elem) == 2:
            raise ValidationError("All elements in Spectrum list must have two items")
        if not all(isinstance(n, (int, float)) for n in elem):
            raise ValidationError("All items in Septrum list elements must be numbers")

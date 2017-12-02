from django.core.exceptions import ValidationError
from proteins.models import Spectrum
import ast


def validate_spectrum(value):
    if not value:
        return None
    if isinstance(value, Spectrum):
        return
    try:
        obj = ast.literal_eval(value)
    except:
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

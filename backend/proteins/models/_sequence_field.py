from django.db import models

from fpseq import FPSeq
from proteins.validators import protein_sequence_validator


class SequenceField(models.CharField):
    def __init__(self, *args, **kwargs):
        kwargs["max_length"] = 1024
        kwargs["validators"] = [protein_sequence_validator]
        super().__init__(*args, **kwargs)

    def deconstruct(self):
        name, path, args, kwargs = super().deconstruct()
        del kwargs["max_length"]
        del kwargs["validators"]
        return name, path, args, kwargs

    def from_db_value(self, value, expression, connection):
        # Skip validation for database values - they're already validated
        return FPSeq(value, validate=False) if value else None

    def to_python(self, value):
        if isinstance(value, FPSeq):
            return value
        # New values should still be validated
        return FPSeq(value, validate=True) if value else None

    def get_prep_value(self, value):
        return str(value) if value else None

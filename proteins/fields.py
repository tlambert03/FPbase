from django.db.models import Lookup, fields, TextField
from django.core.exceptions import ValidationError
import json
from .helpers import wave_to_hex


class Around(Lookup):
    lookup_name = 'around'

    def as_sql(self, compiler, connection):
        lhs, lhs_params = self.process_lhs(compiler, connection)
        rhs, rhs_params = self.process_rhs(compiler, connection)
        params = lhs_params + rhs_params + lhs_params + rhs_params
        return '%s > %s - 16 AND %s < %s + 16' % (lhs, rhs, lhs, rhs), params


fields.PositiveSmallIntegerField.register_lookup(Around)


class NotEqual(Lookup):
    lookup_name = 'ne'

    def as_sql(self, compiler, connection):
        lhs, lhs_params = self.process_lhs(compiler, connection)
        rhs, rhs_params = self.process_rhs(compiler, connection)
        params = lhs_params + rhs_params
        return '%s <> %s' % (lhs, rhs), params


fields.CharField.register_lookup(NotEqual)


class Spectrum(object):
    """ Python class for spectra as a list of lists """

    def __init__(self, data=None):
        if data:
            if not isinstance(data, list):  # must be a list
                raise TypeError("Spectrum object must be of type List")
            if not all(isinstance(elem, list) for elem in data):  # must be list of lists
                raise TypeError("Spectrum object must be a list of lists")
            for elem in data:
                if not len(elem) == 2:
                    raise TypeError("All elements in Spectrum list must have two items")
                if not all(isinstance(n, (int, float)) for n in elem):
                    raise TypeError("All items in Septrum list elements must be numbers")
        self.data = data

    @property
    def x(self):
        self._x = []
        for i in self.data:
            self._x.append(i[0])
        return self._x

    @property
    def y(self):
        self._y = []
        for i in self.data:
            self._y.append(i[1])
        return self._y

    @property
    def color(self):
        return wave_to_hex(self.peak_wave)

    @property
    def peak_wave(self):
        return self.x[self.y.index(max(self.y))]

    @property
    def min_wave(self):
        return self.x[0]

    @property
    def max_wave(self):
        return self.x[-1]

    def __str__(self):
        return json.dumps(self.data)

    def width(self, height=0.5):
        try:
            upindex = next(x[0] for x in enumerate(self.y) if x[1] > height)
            downindex = len(self.y) - next(x[0] for x in enumerate(reversed(self.y)) if x[1] > height)
            return (self.x[upindex], self.x[downindex])
        except Exception:
            return False

    def change_x(self, value):
        if not isinstance(value, list):
            raise TypeError("X values must be a python list")
        if len(value) != len(self.data):
            raise ValueError("Error: array length must match existing data")
        for i in range(len(value)):
            self.data[i][0] = value[i]

    def change_y(self, value):
        if not isinstance(value, list):
            raise TypeError("Y values must be a python list")
        if len(value) != len(self.data):
            raise ValueError("Error: array length must match existing data")
        for i in range(len(value)):
            self.data[i][1] = value[i]

    def nvd3Format(self):
        output = []
        # arrayLength = len(self.data)
        for wave in range(350, int(self.min_wave)):
            output.append({'x': wave, 'y': 0})
        for elem in self.data:
            output.append({'x': elem[0], 'y': elem[1]})
        for wave in range(int(self.max_wave), 751):
            output.append({'x': wave, 'y': 0})
        return output

    def wave_value_pairs(self):
        output = {}
        # arrayLength = len(self.data)
        for elem in self.data:
            output[elem[0]] = elem[1]
        return output


class SpectrumField(TextField):
    description = "Stores a spectrum object"

    def __init__(self, *args, **kwargs):
        super(SpectrumField, self).__init__(*args, **kwargs)

    def from_db_value(self, value, expression, connection, context):
        if not value:
            return None
        return Spectrum(json.loads(value))

    def to_python(self, value):
        if isinstance(value, Spectrum):
            return value

        if not value:
            return None

        try:
            obj = json.loads(value)
            return Spectrum(obj)
        except Exception:
            raise ValidationError("Invalid input for a Spectrum instance")

    def get_prep_value(self, value):
        if value is None:
            return value
        return str(value)


# from django.utils.inspect import func_supports_parameter
# from django.contrib.postgres.fields import ArrayField

# class specField(ArrayField):

#     def _from_db_value(self, value, expression, connection):
#         if value is None:
#             return value
#         return Spectrum([
#             self.base_field.from_db_value(item, expression, connection, {})
#             if func_supports_parameter(self.base_field.from_db_value, 'context')  # RemovedInDjango30Warning
#             else self.base_field.from_db_value(item, expression, connection)
#             for item in value
#         ])

#     def to_python(self, value):
#         if isinstance(value, Spectrum):
#             return value
#         # Assume we're deserializing
#         if isinstance(value, str):
#             # Assume we're deserializing
#             vals = json.loads(value)
#             value = [self.base_field.to_python(val) for val in vals]
#         try:
#             return Spectrum(value)
#         except Exception:
#             raise ValidationError("Invalid input for a Spectrum instance")

#     def get_db_prep_value(self, value, connection, prepared=False):
#         if isinstance(value, Spectrum):
#             value = value.data
#         if isinstance(value, (list, tuple)):
#             return [self.base_field.get_db_prep_value(i, connection, prepared=False) for i in value]
#         return value

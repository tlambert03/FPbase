from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.core.exceptions import ValidationError
from model_utils.models import TimeStampedModel
from .mixins import Authorable

# WORK IN PROGRESS


class SpectrumData(ArrayField):

    def __init__(self, base_field=None, size=None, **kwargs):
        if not base_field:
            base_field = ArrayField(models.FloatField(max_length=10), size=2)
        super().__init__(base_field, size, **kwargs)

    # def from_db_value(self, value, expression, connection, *args, **kwargs):
    #    if value is None:
    #        return value
    #    return value

    # def to_python(self, value):
    #    if isinstance(value, np.ndarray):
    #        return value
    #    return super().to_python(np.array(value))

    # def get_db_prep_value(self, value, connection, prepared=False):
    #    print('prep value')
    #    return super().get_db_prep_value(np.array(value).tolist(), connection, prepared)

    # def value_to_string(self, obj):
    #    print('to str')
    #    return super().value_to_string(obj.tolist())

    def validate(self, value, model_instance):
        super().validate(value, model_instance)
        for elem in value:
            if not len(elem) == 2:
                raise ValidationError("All elements in Spectrum list must have two items")
            if not all(isinstance(n, (int, float)) for n in elem):
                raise ValidationError("All items in Septrum list elements must be numbers")


class Spectrum(Authorable, TimeStampedModel):
    EX = 'ex'
    ABS = 'abs'
    EM = 'em'
    TWOP = 'twop'
    SPECTRA_TYPES = (
        (EX, 'excitation'),
        (ABS, 'absorption'),
        (EM, 'emission'),
        (TWOP, 'two photon absorption'),
    )

    data = SpectrumData()
    stype = models.CharField(max_length=4, choices=SPECTRA_TYPES, blank=True, verbose_name='Spectra Type')
    state = models.ForeignKey('State', related_name="spectra", help_text="The protein(state) to which this spectrum belongs", on_delete=models.CASCADE, blank=True, null=True)

    def __str__(self):
        return "{} {}".format(self.state.protein if self.state else 'unowned', self.stype)

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
    def peak_wave(self):
        return self.x[self.y.index(max(self.y))]

    @property
    def min_wave(self):
        return self.x[0]

    @property
    def max_wave(self):
        return self.x[-1]

    @property
    def step(self):
        s = set()
        for i in range(len(self.x)-1):
            s.add(self.x[i+1] - self.x[i])
        if len(s) > 1:
            return False
        return list(s)[0]

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

    def save(self, *args, **kwargs):
        self.full_clean()
        super().save(*args, **kwargs)

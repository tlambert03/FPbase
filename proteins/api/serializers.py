from ..models import Protein, State, StateTransition
from rest_framework import serializers
from drf_tweaks.serializers import ModelSerializer


class StateTransitionSerializer(serializers.ModelSerializer):
    from_state = serializers.SlugRelatedField(read_only=True, slug_field='slug')
    to_state = serializers.SlugRelatedField(read_only=True, slug_field='slug')

    class Meta:
        model = StateTransition
        fields = ('from_state', 'to_state', 'trans_wave')


class SpectrumField(serializers.Field):
    def to_representation(self, obj):
        return obj.data


class StateSpectraSerializer(serializers.ModelSerializer):
    ex_spectra = SpectrumField()  # adds significant time overhead
    em_spectra = SpectrumField()  # adds significant time overhead

    class Meta:
        model = State
        fields = ('name', 'ex_max', 'ex_spectra', 'em_max', 'em_spectra', 'ext_coeff', 'qy')


class ProteinSpectraSerializer(ModelSerializer):
    states = StateSpectraSerializer(many=True, read_only=True)

    class Meta:
        model = Protein
        fields = ('name', 'slug', 'states')

    def to_representation(self, obj):
        """Move fields from spectra to protein representation."""
        representation = super().to_representation(obj)
        spectra_repr = representation.pop('states')
        representation['spectra'] = []
        for spectrum in spectra_repr:
            if spectrum['ex_spectra']:
                representation['spectra'].append({
                    'state': spectrum['name'] + str('_ex'),
                    'ec': spectrum['ext_coeff'],
                    'max': spectrum['ex_max'],
                    'data': spectrum['ex_spectra'],
                })
            if spectrum['em_spectra']:
                representation['spectra'].append({
                    'state': spectrum['name'] + str('_em'),
                    'qy': spectrum['qy'],
                    'max': spectrum['em_max'],
                    'data': spectrum['em_spectra']
                })
        return representation


class StateSerializer(ModelSerializer):
    protein = serializers.SlugRelatedField(slug_field='slug', read_only=True)

    class Meta:
        model = State
        fields = ('slug', 'protein', 'name', 'ex_max', 'em_max', 'ex_spectra', 'em_spectra',
                  'ext_coeff', 'qy', 'pka', 'maturation', 'lifetime', 'brightness',
                  )
        on_demand_fields = ('protein', 'ex_spectra', 'em_spectra')


class ProteinSerializer(ModelSerializer):
    # url = serializers.CharField(source='get_absolute_url', read_only=True)
    states = StateSerializer(many=True, read_only=True)
    transitions = StateTransitionSerializer(many=True, read_only=True)
    uuid = serializers.UUIDField(format='hex')
    doi = serializers.SlugRelatedField(source='primary_reference', slug_field='doi', read_only=True)

    class Meta:
        model = Protein
        fields = (
            # 'url',
            'uuid',
            'name',
            'slug',
            'seq',
            'ipg_id',
            'genbank',
            'uniprot',
            'pdb',
            'agg',
            'switch_type',
            'states',
            'transitions',
            'doi',
        )
        on_demand_fields = ()


class BleachField(serializers.Field):

    def to_representation(self, obj):
        if obj.default_state.bleach_measurements.count():
            rate = obj.default_state.bleach_measurements.first().rate
            return rate
        else:
            return None


# NOT DRY
# TODO: figure out how to combine this with above
class BasicProteinSerializer(ModelSerializer, serializers.HyperlinkedModelSerializer):
    # states = StateSerializer(many=True, read_only=True)
    url = serializers.CharField(source='get_absolute_url', read_only=True)
    uuid = serializers.UUIDField(format='hex')
    ex_max = serializers.IntegerField(source='default_state.ex_max', read_only=True)
    em_max = serializers.IntegerField(source='default_state.em_max', read_only=True)
    ex_spectra = serializers.CharField(source='default_state.ex_spectra', read_only=True)
    em_spectra = serializers.CharField(source='default_state.em_spectra', read_only=True)
    ext_coeff = serializers.FloatField(source='default_state.ext_coeff', read_only=True)
    qy = serializers.FloatField(source='default_state.qy', read_only=True)
    brightness = serializers.FloatField(source='default_state.brightness', read_only=True)
    bleach = BleachField(source='*', read_only=True)
    maturation = serializers.FloatField(source='default_state.maturation', read_only=True)
    lifetime = serializers.FloatField(source='default_state.lifetime', read_only=True)
    pka = serializers.FloatField(source='default_state.pka', read_only=True)
    stokes = serializers.FloatField(source='default_state.stokes', read_only=True)

    class Meta:
        model = Protein
        fields = (
            'url',
            'uuid',
            'name',
            'stokes',
            'slug',
            'ipg_id',
            'agg',
            'ex_max',
            'ex_spectra',
            'em_max',
            'em_spectra',
            'ext_coeff',
            'qy',
            'pka',
            'brightness',
            'bleach',
            'maturation',
            'lifetime',
        )
        on_demand_fields = ('uuid', 'ex_spectra', 'em_spectra')

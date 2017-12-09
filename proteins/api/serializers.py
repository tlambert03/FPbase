from ..models import Protein, State, StateTransition
from rest_framework import serializers


class StateTransitionSerializer(serializers.ModelSerializer):
    from_state = serializers.SlugRelatedField(read_only=True, slug_field='slug')
    to_state = serializers.SlugRelatedField(read_only=True, slug_field='slug')

    class Meta:
        model = StateTransition
        fields = ('from_state', 'to_state', 'trans_wave')


class StateSerializer(serializers.ModelSerializer):
    class Meta:
        model = State
        fields = ('slug', 'ex_max', 'em_max', 'ex_spectra', 'em_spectra',
                  'ext_coeff', 'qy', 'pka', 'maturation', 'lifetime'
                  )


class ProteinSerializer(serializers.HyperlinkedModelSerializer):
    # url = serializers.CharField(source='get_absolute_url', read_only=True)
    states = StateSerializer(many=True, read_only=True)
    transitions = StateTransitionSerializer(many=True, read_only=True)
    uuid = serializers.UUIDField(format='hex')

    class Meta:
        model = Protein
        fields = (
            # 'url',
            'uuid',
            'name',
            'slug',
            'seq',
            'gb_prot',
            'ipg_id',
            'agg',
            'switch_type',
            'states',
            'transitions'
        )


# NOT DRY
# TODO: figure out how to combine this with above
class BasicProteinSerializer(serializers.HyperlinkedModelSerializer):
    # states = StateSerializer(many=True, read_only=True)
    url = serializers.CharField(source='get_absolute_url', read_only=True)
    lambda_ex = serializers.IntegerField(source='default_state.ex_max', read_only=True)
    lambda_em = serializers.IntegerField(source='default_state.em_max', read_only=True)
    E = serializers.FloatField(source='default_state.ext_coeff', read_only=True)
    QY = serializers.FloatField(source='default_state.qy', read_only=True)
    brightness = serializers.FloatField(source='default_state.brightness', read_only=True)
    bleach = serializers.FloatField(source='default_state.bleach', read_only=True)
    mature = serializers.FloatField(source='default_state.maturation', read_only=True)
    lifetime = serializers.FloatField(source='default_state.lifetime', read_only=True)
    pka = serializers.FloatField(source='default_state.pka', read_only=True)
    stokes = serializers.FloatField(source='default_state.stokes', read_only=True)

    class Meta:
        model = Protein
        fields = (
            'url',
            'name',
            'stokes',
            'slug',
            'ipg_id',
            'agg',
            'lambda_ex',
            'lambda_em',
            'E',
            'QY',
            'pka',
            'brightness',
            'bleach',
            'mature',
            'lifetime',
        )

from rest_framework import serializers

from ..models import Protein, Spectrum, State, StateTransition
from ._tweaks import ModelSerializer


class SpectrumSerializer(serializers.ModelSerializer):
    owner_id = serializers.IntegerField(source="owner.id")
    owner_slug = serializers.CharField(source="owner.slug")
    protein_name = serializers.SerializerMethodField()
    protein_slug = serializers.SerializerMethodField()

    class Meta:
        model = Spectrum
        fields = (
            "data",
            "category",
            "subtype",
            "id",
            "name",
            "owner_id",
            "owner_slug",
            "protein_name",
            "protein_slug",
            "color",
            "min_wave",
            "max_wave",
            "peak_wave",
        )

    def get_protein_name(self, obj):
        if obj.owner_state:
            return obj.owner_state.protein.name

    def get_protein_slug(self, obj):
        if obj.owner_state:
            return obj.owner_state.protein.slug


class StateTransitionSerializer(serializers.ModelSerializer):
    from_state = serializers.SlugRelatedField(read_only=True, slug_field="slug")
    to_state = serializers.SlugRelatedField(read_only=True, slug_field="slug")

    class Meta:
        model = StateTransition
        fields = ("from_state", "to_state", "trans_wave")


class SpectrumField(serializers.Field):
    def to_representation(self, obj):
        return obj.data


class StateSpectraSerializer(serializers.ModelSerializer):
    ex_spectrum = SpectrumField()  # adds significant time overhead
    em_spectrum = SpectrumField()  # adds significant time overhead

    class Meta:
        model = State
        fields = (
            "name",
            "ex_max",
            "ex_spectrum",
            "em_max",
            "em_spectrum",
            "ext_coeff",
            "qy",
        )


class ProteinSpectraSerializer(ModelSerializer):
    states = StateSpectraSerializer(many=True, read_only=True)

    class Meta:
        model = Protein
        fields = ("name", "slug", "states")

    def to_representation(self, obj):
        """Move fields from spectra to protein representation."""
        representation = super().to_representation(obj)
        spectra_repr = representation.pop("states")
        representation["spectra"] = []
        for spectrum in spectra_repr:
            if spectrum["ex_spectrum"]:
                representation["spectra"].append(
                    {
                        "state": spectrum["name"] + "_ex",
                        "ec": spectrum["ext_coeff"],
                        "max": spectrum["ex_max"],
                        "data": spectrum["ex_spectrum"],
                    }
                )
            if spectrum["em_spectrum"]:
                representation["spectra"].append(
                    {
                        "state": spectrum["name"] + "_em",
                        "qy": spectrum["qy"],
                        "max": spectrum["em_max"],
                        "data": spectrum["em_spectrum"],
                    }
                )
        return representation


class StateSerializer(ModelSerializer):
    protein = serializers.SlugRelatedField(slug_field="slug", read_only=True)

    class Meta:
        model = State
        fields = (
            "slug",
            "protein",
            "name",
            "ex_max",
            "em_max",
            "ex_spectrum",
            "em_spectrum",
            "ext_coeff",
            "qy",
            "pka",
            "maturation",
            "lifetime",
            "brightness",
        )
        on_demand_fields = ("protein", "ex_spectrum", "em_spectrum")


class ProteinSerializer(ModelSerializer):
    # url = serializers.CharField(source='get_absolute_url', read_only=True)
    states = StateSerializer(many=True, read_only=True)
    transitions = StateTransitionSerializer(many=True, read_only=True)
    doi = serializers.SlugRelatedField(source="primary_reference", slug_field="doi", read_only=True)

    class Meta:
        model = Protein
        fields = (
            # 'url',
            "uuid",
            "name",
            "slug",
            "seq",
            "ipg_id",
            "genbank",
            "uniprot",
            "pdb",
            "agg",
            "switch_type",
            "states",
            "transitions",
            "doi",
        )
        on_demand_fields = ()


class ProteinSerializer2(ModelSerializer):
    states = serializers.SlugRelatedField(many=True, read_only=True, slug_field="slug")
    transitions = serializers.IntegerField(source="transitions.count", read_only=True)
    doi = serializers.SlugRelatedField(source="primary_reference", slug_field="doi", read_only=True)

    class Meta:
        model = Protein
        fields = (
            "name",
            "slug",
            "seq",
            "agg",
            "doi",
            "states",
            "pdb",
            "switch_type",
            "genbank",
            "uniprot",
            "ipg_id",
            "transitions",
        )
        on_demand_fields = (
            "pdb",
            "switch_type",
            "genbank",
            "uniprot",
            "ipg_id",
            "seq",
            "transitions",
        )


# NOT DRY
# TODO: figure out how to combine this with above
class BasicProteinSerializer(ModelSerializer, serializers.HyperlinkedModelSerializer):
    # states = StateSerializer(many=True, read_only=True)
    url = serializers.CharField(source="get_absolute_url", read_only=True)
    ex_max = serializers.IntegerField(source="default_state.ex_max", read_only=True)
    em_max = serializers.IntegerField(source="default_state.em_max", read_only=True)
    ex_spectrum = serializers.CharField(source="default_state.ex_spectrum", read_only=True)
    em_spectrum = serializers.CharField(source="default_state.em_spectrum", read_only=True)
    ext_coeff = serializers.FloatField(source="default_state.ext_coeff", read_only=True)
    qy = serializers.FloatField(source="default_state.qy", read_only=True)
    brightness = serializers.FloatField(source="default_state.brightness", read_only=True)
    bleach = serializers.FloatField(source="rate", read_only=True)
    maturation = serializers.FloatField(source="default_state.maturation", read_only=True)
    lifetime = serializers.FloatField(source="default_state.lifetime", read_only=True)
    pka = serializers.FloatField(source="default_state.pka", read_only=True)
    stokes = serializers.FloatField(source="default_state.stokes", read_only=True)

    class Meta:
        model = Protein
        fields = (
            "url",
            "uuid",
            "name",
            "stokes",
            "slug",
            "ipg_id",
            "agg",
            "ex_max",
            "ex_spectrum",
            "em_max",
            "em_spectrum",
            "ext_coeff",
            "qy",
            "pka",
            "brightness",
            "bleach",
            "maturation",
            "lifetime",
            "cofactor",
        )
        on_demand_fields = ("uuid", "ex_spectrum", "em_spectrum")

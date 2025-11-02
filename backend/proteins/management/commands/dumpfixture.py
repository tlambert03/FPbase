# proteins/management/commands/dump_proteins_with_relations.py

import itertools
import math
from collections.abc import Iterable, Sequence

from django.core.management.base import BaseCommand
from django.core.serializers.json import Serializer as JSONSerializer

from proteins.models import FilterPlacement, Microscope, OpticalConfig, Protein

DEST = "backend/tests_e2e/fixtures/test_data.json"


class CleanJSONSerializer(JSONSerializer):
    def get_dump_object(self, obj):
        data = super().get_dump_object(obj)
        for field in ("email", "password", "groups"):
            data["fields"].pop(field, None)
        data["fields"] = self._clean_nan(data["fields"])
        return data

    def _clean_nan(self, obj):
        """Recursively replace NaN with None"""
        if isinstance(obj, float):
            if math.isnan(obj) or math.isinf(obj):
                return None
            return obj
        elif isinstance(obj, dict):
            return {k: self._clean_nan(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [self._clean_nan(item) for item in obj]
        return obj


def _gather_protein(proteins: Sequence[Protein]) -> Iterable[object]:
    for protein in proteins:
        yield protein
        yield protein.created_by

        for state in protein.states.all():
            yield state
            yield from state.spectra.all()
        yield from protein.references.all()
        yield protein.primary_reference
        yield protein.parent_organism


def _gather_scope(scopes: Sequence[Microscope]) -> Iterable[object]:
    for scope in scopes:
        yield scope
        yield from scope.extra_lights.all()
        yield from scope.extra_cameras.all()
        for oc in OpticalConfig.objects.filter(microscope=scope):
            yield oc
            if oc.light:
                yield oc.light
            if oc.camera:
                yield oc.camera

            for fp in FilterPlacement.objects.filter(config=oc):
                yield fp
                yield fp.filter


class Command(BaseCommand):
    help = "Dump proteins with all related objects"

    def add_arguments(self, parser):
        parser.add_argument("--indent", type=int, default=2)
        parser.add_argument("--dest", type=str, default=DEST)

    def handle(self, *args, **options):
        STANDARD_PROTS = ["EGFP", "mCherry"]
        proteins = Protein.objects.filter(name__in=STANDARD_PROTS).prefetch_related(
            "states__spectra",
            "references",
            "parent_organism",
        )
        scopes = Microscope.objects.filter(
            pk__in=[
                "wKqWbgApvguSNDSRZNSfpN",  # Example Simple Widefield
                "4yL4ggAozzcMwTU4Ae7zxF",  # Example Yokogawa Setup
            ]
        )

        # Remove duplicates while preserving order
        seen = set()
        unique_objects = []
        for obj in itertools.chain(_gather_protein(proteins), _gather_scope(scopes)):
            key = (type(obj), obj.pk)
            if key not in seen:
                seen.add(key)
                unique_objects.append(obj)

        # Serialize
        data = CleanJSONSerializer().serialize(
            unique_objects,
            indent=options["indent"],
            stream=None,
            fields=None,
        )

        with open(options["dest"], "w") as f:
            f.write(data)

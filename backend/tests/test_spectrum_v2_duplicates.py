"""Test duplicate spectrum validation."""

import json

import pytest
from django import forms

from proteins.forms.spectrum_v2 import _validate_spectrum_json


@pytest.mark.django_db
class TestDuplicateSpectraValidation:
    """Test that duplicate spectra within the same form are rejected."""

    def test_duplicate_spectra_rejected(self):
        """Test that duplicate spectra (same category/owner/subtype) are rejected."""
        spectra_json = json.dumps(
            [
                {
                    "data": [[400, 0.1], [500, 1.0], [600, 0.5]],
                    "category": "d",
                    "owner": "Alexa Fluor 488",
                    "subtype": "ex",
                    "scale_factor": None,
                    "ph": None,
                    "solvent": None,
                    "peak_wave": 500,
                    "column_name": "Column A",
                },
                {
                    "data": [[400, 0.2], [500, 0.9], [600, 0.4]],  # Different data
                    "category": "d",
                    "owner": "Alexa Fluor 488",  # Same owner
                    "subtype": "ex",  # Same subtype
                    "scale_factor": None,
                    "ph": None,
                    "solvent": None,
                    "peak_wave": 500,
                    "column_name": "Column B",  # Different column name
                },
            ]
        )

        with pytest.raises(forms.ValidationError, match="Duplicate spectrum detected"):
            _validate_spectrum_json(spectra_json)

    def test_different_subtypes_allowed(self):
        """Test that same owner with different subtypes is allowed."""
        spectra_json = json.dumps(
            [
                {
                    "data": [[400, 0.1], [500, 1.0], [600, 0.5]],
                    "category": "d",
                    "owner": "Alexa Fluor 488",
                    "subtype": "ex",
                    "scale_factor": None,
                    "ph": None,
                    "solvent": None,
                    "peak_wave": 500,
                    "column_name": "Excitation",
                },
                {
                    "data": [[500, 0.2], [520, 0.9], [600, 0.4]],
                    "category": "d",
                    "owner": "Alexa Fluor 488",  # Same owner
                    "subtype": "em",  # Different subtype - should be OK
                    "scale_factor": None,
                    "ph": None,
                    "solvent": None,
                    "peak_wave": 520,
                    "column_name": "Emission",
                },
            ]
        )

        # Should not raise
        result = _validate_spectrum_json(spectra_json)
        assert len(result) == 2

    def test_case_insensitive_duplicate_detection(self):
        """Test that duplicate detection is case-insensitive."""
        spectra_json = json.dumps(
            [
                {
                    "data": [[400, 0.1], [500, 1.0], [600, 0.5]],
                    "category": "d",
                    "owner": "Alexa Fluor 488",
                    "subtype": "ex",
                    "scale_factor": None,
                    "ph": None,
                    "solvent": None,
                    "peak_wave": 500,
                    "column_name": "Column A",
                },
                {
                    "data": [[400, 0.2], [500, 0.9], [600, 0.4]],
                    "category": "d",
                    "owner": "ALEXA FLUOR 488",  # Different case
                    "subtype": "ex",
                    "scale_factor": None,
                    "ph": None,
                    "solvent": None,
                    "peak_wave": 500,
                    "column_name": "Column B",
                },
            ]
        )

        with pytest.raises(forms.ValidationError, match="Duplicate spectrum detected"):
            _validate_spectrum_json(spectra_json)

from __future__ import annotations

import pytest
from django.core.exceptions import ValidationError

from proteins.models import Spectrum


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_spectrum_with_non_finite_values_is_invalid(bad: float) -> None:
    spectrum = Spectrum(category=Spectrum.FILTER, subtype=Spectrum.BS)
    spectrum.data = [[400, 0.5], [401, bad], [402, 0.7]]
    with pytest.raises(ValidationError, match="finite"):
        spectrum.clean()


def test_spectrum_with_finite_values_is_valid() -> None:
    spectrum = Spectrum(category=Spectrum.FILTER, subtype=Spectrum.BS)
    spectrum.data = [[400, 0.5], [401, 0.6], [402, 0.7]]
    spectrum.clean()

from __future__ import annotations

import pytest

from proteins.models import Spectrum
from proteins.util.spectra import interp_linear


@pytest.mark.parametrize(
    ("first", "last", "expected"),
    [
        (300.0, 310.0, range(300, 311)),  # ends on a whole nm: keep it
        (300.5, 309.5, range(301, 310)),
    ],
)
def test_interp_linear_covers_every_whole_nm(first: float, last: float, expected: range):
    waves = [first + i / 2 for i in range(int((last - first) * 2) + 1)]
    values = [w / 1000 for w in waves]
    x, y = interp_linear(waves, values)
    assert list(x) == list(expected)
    assert y == pytest.approx([w / 1000 for w in expected])


def test_resampled_spectrum_keeps_its_last_wavelength():
    waves = [400 + i / 2 for i in range(201)]  # 400 to 500 in 0.5 nm steps
    spectrum = Spectrum()
    spectrum._set_spectrum_data(waves, [0.5] * len(waves))
    assert (spectrum.min_wave, spectrum.max_wave) == (400, 500)
    assert len(spectrum.y) == 101

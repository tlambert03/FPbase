from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
from django.core.exceptions import ValidationError
from django.db import models
from model_utils.models import TimeStampedModel

from proteins.models.mixins import Authorable
from references.models import Reference

if TYPE_CHECKING:
    from proteins.models import BleachMeasurement, State


class CurveData(Authorable, TimeStampedModel):
    """Abstract base for 2D curve data with arbitrary X-axis spacing.

    Unlike Spectrum (which assumes 1nm X-axis spacing), this model stores
    both X and Y values explicitly to support non-uniform sampling.

    Subclasses define the curve type and ownership relationship.
    """

    # Binary storage for X and Y values (float32)
    x_values = models.BinaryField(help_text="X values as float32 binary array")
    y_values = models.BinaryField(help_text="Y values as float32 binary array")

    # Cached bounds (computed on save)
    x_min = models.FloatField(null=True, blank=True, help_text="Minimum X value")
    x_max = models.FloatField(null=True, blank=True, help_text="Maximum X value")

    # Provenance
    reference_id: int | None
    reference: models.ForeignKey[Reference | None] = models.ForeignKey(
        Reference,
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="%(class)s_curves",
    )
    source = models.CharField(max_length=128, blank=True, help_text="Source of the curve data")
    notes = models.TextField(blank=True, help_text="Additional notes or conditions")

    class Meta:
        abstract = True

    def __str__(self) -> str:
        return f"{self.__class__.__name__} ({self.x_min:.2f}-{self.x_max:.2f})"

    def save(self, *args, **kwargs) -> None:
        # Cache min/max from data
        if self.x_values:
            x_arr = self._decode_values(self.x_values)
            self.x_min = float(x_arr.min())
            self.x_max = float(x_arr.max())
        super().save(*args, **kwargs)

    def clean(self) -> None:
        super().clean()
        if self.x_values and self.y_values:
            x_len = len(self.x_values) // 4  # float32 = 4 bytes
            y_len = len(self.y_values) // 4
            if x_len != y_len:
                raise ValidationError(
                    f"X values length ({x_len}) must match Y values length ({y_len})"
                )

    # =========================================================================
    # Binary encoding/decoding (following Spectrum pattern)
    # =========================================================================

    @staticmethod
    def _encode_values(values: list[float]) -> bytes:
        """Encode list of floats to float32 binary data."""
        return np.asarray(values, dtype="<f4").tobytes()

    @staticmethod
    def _decode_values(data: bytes | memoryview) -> np.ndarray:
        """Decode float32 binary data to numpy array."""
        if isinstance(data, memoryview):
            data = bytes(data)
        return np.frombuffer(data, dtype="<f4")

    @cached_property
    def x(self) -> list[float]:
        """X values decoded from binary storage."""
        if not self.x_values:
            return []
        return self._decode_values(self.x_values).tolist()

    @cached_property
    def y(self) -> list[float]:
        """Y values decoded from binary storage."""
        if not self.y_values:
            return []
        return self._decode_values(self.y_values).tolist()

    def set_data(self, x_values: list[float], y_values: list[float]) -> None:
        """Set curve data from X and Y value lists.

        Parameters
        ----------
        x_values
            X-axis values (pH for PKA, time for bleach/maturation).
        y_values
            Y-axis values (typically normalized intensity 0-1).
        """
        if len(x_values) != len(y_values):
            raise ValueError("X and Y arrays must have same length")

        # Clear cached properties
        for prop in ("x", "y"):
            if prop in self.__dict__:
                del self.__dict__[prop]

        self.x_values = self._encode_values(x_values)
        self.y_values = self._encode_values(y_values)

    @property
    def data(self) -> list[tuple[float, float]]:
        """Return data as list of (x, y) tuples."""
        return list(zip(self.x, self.y))

    @data.setter
    def data(self, value: list[list[float]] | list[tuple[float, float]]) -> None:
        """Set data from list of [x, y] pairs."""
        if not value:
            return
        x_vals = [point[0] for point in value]
        y_vals = [point[1] for point in value]
        self.set_data(x_vals, y_vals)

    def d3data(self) -> list[dict[str, float]]:
        """Return curve data in D3.js format."""
        return [{"x": x, "y": y} for x, y in zip(self.x, self.y)]


class PKACurve(CurveData):
    """pH sensitivity curve (intensity vs pH).

    X-axis: pH values (typically 4-10)
    Y-axis: normalized fluorescence intensity (0-1)
    """

    state_id: int
    state: models.ForeignKey[State] = models.ForeignKey(
        "State",
        on_delete=models.CASCADE,
        related_name="pka_curves",
        help_text="Protein state this curve belongs to",
    )

    class Meta:
        verbose_name = "pKa Curve"
        verbose_name_plural = "pKa Curves"
        indexes = [
            models.Index(fields=["state_id", "status"], name="pkacurve_state_status_idx"),
        ]

    def __str__(self) -> str:
        return f"{self.state} pKa curve"


class BleachCurve(CurveData):
    """Photobleaching decay curve (intensity vs time).

    X-axis: time in seconds
    Y-axis: normalized fluorescence intensity (0-1)

    Links to BleachMeasurement which contains experimental conditions
    (power, modality, temperature, etc.).
    """

    bleach_measurement_id: int
    bleach_measurement: models.ForeignKey[BleachMeasurement] = models.ForeignKey(
        "BleachMeasurement",
        on_delete=models.CASCADE,
        related_name="curves",
        help_text="BleachMeasurement this curve belongs to",
    )

    class Meta:
        verbose_name = "Bleach Curve"
        verbose_name_plural = "Bleach Curves"
        indexes = [
            models.Index(
                fields=["bleach_measurement_id", "status"],
                name="bleachcurve_meas_status_idx",
            ),
        ]

    def __str__(self) -> str:
        return f"{self.bleach_measurement.state} bleach curve"


class MaturationCurve(CurveData):
    """Maturation time course curve (intensity vs time).

    X-axis: time in minutes
    Y-axis: normalized fluorescence intensity (0-1)
    """

    state_id: int
    state: models.ForeignKey[State] = models.ForeignKey(
        "State",
        on_delete=models.CASCADE,
        related_name="maturation_curves",
        help_text="Protein state this curve belongs to",
    )
    temperature = models.FloatField(
        null=True,
        blank=True,
        help_text="Temperature in °C during maturation measurement",
    )

    class Meta:
        verbose_name = "Maturation Curve"
        verbose_name_plural = "Maturation Curves"
        indexes = [
            models.Index(fields=["state_id", "status"], name="matcurve_state_status_idx"),
        ]

    def __str__(self) -> str:
        return f"{self.state} maturation curve"

"""Tests for Celery tasks in proteins app."""

from unittest.mock import MagicMock, patch

import pytest

from proteins.tasks import calc_fret, calculate_scope_report


@pytest.mark.django_db
class TestCalculateScopeReport:
    """Test calculate_scope_report task memory optimization."""

    def test_calculate_scope_report_with_outdated_ids(self):
        """Test that calculate_scope_report processes outdated IDs in batches."""
        # Create a list of mock outdated IDs
        outdated_ids = list(range(150))  # 150 IDs to test batching

        with (
            patch("proteins.models.OcFluorEff") as mock_eff_model,
            patch("proteins.models.Microscope") as mock_microscope,
        ):
            # Mock the filter().filter() chain
            mock_queryset = MagicMock()
            mock_queryset.__iter__ = MagicMock(return_value=iter([]))
            mock_eff_model.objects.filter.return_value = mock_queryset

            # Mock microscope
            mock_microscope.objects.get.return_value = MagicMock()

            # Call the task with .delay() to test batching logic
            # Since we're using mocks, it won't actually execute async
            with patch.object(calculate_scope_report, "update_state"):
                calculate_scope_report.run(scope_id=1, outdated_ids=outdated_ids)

            # Verify that batching occurred (150 IDs / 50 per batch = 3 batches)
            assert mock_eff_model.objects.filter.call_count == 3

    def test_calculate_scope_report_uses_values_list(self):
        """Test that calculate_scope_report uses values_list for memory efficiency."""
        with (
            patch("proteins.models.State") as mock_state,
            patch("proteins.models.Dye") as mock_dye,
            patch("proteins.models.Microscope") as mock_microscope,
        ):
            # Mock the with_spectra().values_list() chain
            mock_state_qs = MagicMock()
            mock_state_qs.values_list.return_value = []
            mock_state.objects.with_spectra.return_value = mock_state_qs

            mock_dye_qs = MagicMock()
            mock_dye_qs.values_list.return_value = []
            mock_dye.objects.with_spectra.return_value = mock_dye_qs

            # Mock microscope
            mock_microscope_instance = MagicMock()
            mock_microscope_instance.optical_configs.count.return_value = 0
            mock_microscope.objects.get.return_value = mock_microscope_instance

            # Call the task directly via .run()
            calculate_scope_report.run(scope_id=1)

            # Verify values_list was called with id and flat=True (memory optimization)
            mock_state_qs.values_list.assert_called_once_with("id", flat=True)
            mock_dye_qs.values_list.assert_called_once_with("id", flat=True)

    def test_calculate_scope_report_signature(self):
        """Test that calculate_scope_report accepts expected parameters."""
        # This is a smoke test to ensure the function signature hasn't changed
        with (
            patch("proteins.models.State") as mock_state,
            patch("proteins.models.Dye") as mock_dye,
            patch("proteins.models.Microscope") as mock_microscope,
            patch("proteins.models.OcFluorEff") as mock_eff_model,
        ):
            mock_state.objects.with_spectra.return_value.values_list.return_value = []
            mock_dye.objects.with_spectra.return_value.values_list.return_value = []
            mock_microscope_instance = MagicMock()
            mock_microscope_instance.optical_configs.count.return_value = 0
            mock_microscope.objects.get.return_value = mock_microscope_instance
            mock_eff_model.objects.filter.return_value = MagicMock(__iter__=lambda self: iter([]))

            # Should not raise any exceptions
            calculate_scope_report.run(scope_id=1)
            calculate_scope_report.run(scope_id=1, outdated_ids=[1, 2, 3])
            calculate_scope_report.run(scope_id=1, fluor_collection=None)


class TestCalcFret:
    """Test calc_fret task."""

    def test_calc_fret_calls_forster_list(self):
        """Test that calc_fret calls forster_list function."""
        with patch("proteins.tasks.forster_list") as mock_forster_list:
            mock_forster_list.return_value = []

            result = calc_fret()

            mock_forster_list.assert_called_once()
            assert result == []

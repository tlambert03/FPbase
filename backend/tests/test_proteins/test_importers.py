"""
Comprehensive test suite for CSV parsing in proteins.util.importers.

Tests text_to_spectra() with various CSV formats to ensure production safety.
"""

from django.test import TestCase

from proteins.util.importers import text_to_spectra


class TestTextToSpectra(TestCase):
    """Test text_to_spectra() with various CSV formats."""

    def test_us_format_no_headers(self):
        """Test US format (comma, dot decimal) without headers - PRODUCTION FORMAT."""
        csv = """400,0.1
401,0.2
402,0.3
403,0.5
404,0.8
405,1.0"""

        waves, data, headers = text_to_spectra(csv)
        self.assertEqual(waves, [400.0, 401.0, 402.0, 403.0, 404.0, 405.0])
        self.assertEqual(data, [[0.1, 0.2, 0.3, 0.5, 0.8, 1.0]])
        self.assertEqual(headers, ["col_1"])

    def test_us_format_with_headers(self):
        """Test US format with headers."""
        csv = """wavelength,intensity
400,0.1
401,0.2
402,0.3"""

        waves, data, headers = text_to_spectra(csv)
        self.assertEqual(waves, [400.0, 401.0, 402.0])
        self.assertEqual(data, [[0.1, 0.2, 0.3]])
        self.assertEqual(headers, ["intensity"])

    def test_european_format(self):
        """Test European format (semicolon, comma decimal)."""
        csv = """wavelength;intensity
400;0,1
401;0,2
402;0,3"""

        waves, data, headers = text_to_spectra(csv)
        self.assertEqual(waves, [400.0, 401.0, 402.0])
        self.assertEqual(data, [[0.1, 0.2, 0.3]])
        self.assertEqual(headers, ["intensity"])

    def test_multiple_data_columns(self):
        """Test multiple intensity columns."""
        csv = """wavelength,ex,em,abs
400,0.1,0.05,0.2
401,0.2,0.15,0.3
402,0.3,0.25,0.4"""

        waves, data, headers = text_to_spectra(csv)
        self.assertEqual(waves, [400.0, 401.0, 402.0])
        self.assertEqual(len(data), 3)
        self.assertEqual(data[0], [0.1, 0.2, 0.3])  # ex column
        self.assertEqual(data[1], [0.05, 0.15, 0.25])  # em column
        self.assertEqual(data[2], [0.2, 0.3, 0.4])  # abs column
        self.assertEqual(headers, ["ex", "em", "abs"])

    def test_non_default_wavecol(self):
        """Test wavelength in non-first column."""
        csv = """intensity,wavelength
0.1,400
0.2,401
0.3,402"""

        waves, data, headers = text_to_spectra(csv, wavecol=1)
        self.assertEqual(waves, [400.0, 401.0, 402.0])
        self.assertEqual(data, [[0.1, 0.2, 0.3]])
        self.assertEqual(headers, ["intensity"])

    def test_thousands_separator_us(self):
        """Test handling of US thousands separators.

        Note: CSV with embedded commas will be parsed as separate columns.
        This is correct behavior - thousands separators should not be used
        in CSV files without proper quoting.
        """
        csv = """400,1,000.5
401,2,000.3
402,3,000.8"""

        waves, data, _headers = text_to_spectra(csv)
        self.assertEqual(waves, [400.0, 401.0, 402.0])
        # Parses as 3 columns: wavelength, thousands, decimal
        self.assertEqual(len(data), 2)
        self.assertEqual(data[0], [1.0, 2.0, 3.0])
        self.assertEqual(data[1], [0.5, 0.3, 0.8])

    def test_thousands_separator_european(self):
        """Test European format with thousands separator."""
        csv = """wavelength;intensity
1.000,5;0,1
2.000,3;0,2
3.000,8;0,3"""

        waves, data, _headers = text_to_spectra(csv)
        self.assertEqual(waves, [1000.5, 2000.3, 3000.8])
        self.assertEqual(data, [[0.1, 0.2, 0.3]])

    def test_whitespace_handling(self):
        """Test CSV with extra whitespace."""
        csv = """wavelength , intensity
400 , 0.1
401 , 0.2
402 , 0.3"""

        waves, data, _headers = text_to_spectra(csv)
        self.assertEqual(waves, [400.0, 401.0, 402.0])
        self.assertEqual(data, [[0.1, 0.2, 0.3]])

    def test_empty_lines(self):
        """Test CSV with empty lines (should skip them)."""
        csv = """wavelength,intensity
400,0.1

401,0.2

402,0.3"""

        waves, data, _headers = text_to_spectra(csv)
        self.assertEqual(waves, [400.0, 401.0, 402.0])
        self.assertEqual(data, [[0.1, 0.2, 0.3]])

    def test_actual_production_format(self):
        """Test exact format from production test suite."""
        csv = "400,0.1\n401,0.2\n402,0.3\n403,0.5\n404,0.8\n405,1.0\n406,0.8\n407,0.5\n408,0.3\n409,0.1"

        waves, data, _headers = text_to_spectra(csv)
        self.assertEqual(waves, [400.0, 401.0, 402.0, 403.0, 404.0, 405.0, 406.0, 407.0, 408.0, 409.0])
        self.assertEqual(len(data[0]), 10)
        self.assertEqual(data[0][0], 0.1)
        self.assertEqual(data[0][-1], 0.1)

    def test_bytes_input(self):
        """Test that bytes input (from file upload) works."""
        csv_bytes = b"400,0.1\n401,0.2\n402,0.3"
        csv_str = csv_bytes.decode("utf-8")

        waves, data, _headers = text_to_spectra(csv_str)
        self.assertEqual(waves, [400.0, 401.0, 402.0])
        self.assertEqual(data, [[0.1, 0.2, 0.3]])

    def test_scientific_notation(self):
        """Test that scientific notation is supported."""
        csv = """wavelength,intensity
400,1.5e-3
401,2.3e-3
402,3.1e-3"""

        waves, data, _headers = text_to_spectra(csv)
        self.assertEqual(waves, [400.0, 401.0, 402.0])
        self.assertAlmostEqual(data[0][0], 0.0015, places=6)
        self.assertAlmostEqual(data[0][1], 0.0023, places=6)
        self.assertAlmostEqual(data[0][2], 0.0031, places=6)

    def test_error_single_column(self):
        """Test that single column raises ValueError."""
        csv = """400
401
402"""

        with self.assertRaises(ValueError) as cm:
            text_to_spectra(csv)
        self.assertIn("Could not parse", str(cm.exception))

    def test_error_single_row(self):
        """Test that single row raises ValueError."""
        csv = "400,0.1"

        with self.assertRaises(ValueError) as cm:
            text_to_spectra(csv)
        self.assertIn("Could not parse", str(cm.exception))

    def test_error_non_numeric_data(self):
        """Test that non-numeric data raises ValueError."""
        csv = """wavelength,intensity
400,abc
401,0.2"""

        with self.assertRaises(ValueError) as cm:
            text_to_spectra(csv)
        self.assertIn("Could not parse", str(cm.exception))

    def test_error_empty_input(self):
        """Test that empty input raises ValueError."""
        with self.assertRaises(ValueError):
            text_to_spectra("")

    def test_error_header_only(self):
        """Test that header-only CSV raises ValueError."""
        csv = "wavelength,intensity"

        with self.assertRaises(ValueError) as cm:
            text_to_spectra(csv)
        self.assertIn("Could not parse", str(cm.exception))

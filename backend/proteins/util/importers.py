from io import StringIO


def read_csv_text(text):
    """Read CSV text and parse to list of numeric rows.

    Tries European format (semicolon separator, comma decimal) first,
    then falls back to US format (comma separator, dot decimal).

    Returns:
        tuple: (headers, data) where headers is list of column names
               and data is list of lists of floats
    """
    import csv

    # Try European format first (semicolon separator, comma decimal)
    # Then fall back to US format (comma separator, dot decimal)
    for delimiter, decimal_char in [(";", ","), (",", ".")]:
        reader = csv.reader(StringIO(text), delimiter=delimiter)
        rows = list(reader)

        if len(rows) < 2:  # Need at least header + 1 data row
            continue

        # Check if first row looks like headers (has non-numeric values)
        first_row = rows[0]
        has_headers = False
        try:
            # Try to convert first row to floats
            for cell in first_row:
                if not cell.strip():
                    continue
                normalized = cell.strip()
                if decimal_char == ",":
                    normalized = normalized.replace(".", "").replace(",", ".")
                else:
                    normalized = normalized.replace(",", "")
                float(normalized)
        except (ValueError, AttributeError):
            has_headers = True

        if has_headers:
            headers = first_row
            data_rows = rows[1:]
        else:
            headers = [f"col_{i}" for i in range(len(first_row))]
            data_rows = rows

        # Try to parse all data rows as floats
        try:
            parsed_data = []
            for row in data_rows:
                if not row or all(not cell.strip() for cell in row):
                    continue
                parsed_row = []
                for cell in row:
                    if not cell.strip():
                        continue
                    # Handle thousands and decimal separators
                    normalized = cell.strip()
                    if decimal_char == ",":
                        normalized = normalized.replace(".", "")  # Remove thousands
                        normalized = normalized.replace(",", ".")  # Decimal comma to dot
                    else:
                        normalized = normalized.replace(",", "")  # Remove thousands
                    parsed_row.append(float(normalized))
                parsed_data.append(parsed_row)

            # Validate: all rows have same length and at least 2 columns
            if not parsed_data:
                continue
            row_lengths = [len(row) for row in parsed_data]
            if len(set(row_lengths)) != 1 or row_lengths[0] < 2:
                continue

            return headers, parsed_data

        except (ValueError, IndexError):
            continue

    raise ValueError("Could not parse text as valid spectra csv.")


def text_to_spectra(text, wavecol=0):
    """Convert text string into spectral data

    Args:
        text (str): string containing csv (usually) data
        wavecol (:obj:`int`, optional): column in csv that contains wavelengths

    Returns:
        tuple: (waves, outdata, headers).  waves is 1D, outdata is MxN, where M
            is the number of data columns and N is the number of wavelengths.
            headers is 1D of length M, containing titles of data columns
    """
    headers, data = read_csv_text(text)

    # Extract wavelength column
    waves = [row[wavecol] for row in data]

    # Extract headers (excluding wavecol)
    out_headers = [str(h) for i, h in enumerate(headers) if i != wavecol]

    # Extract data columns (excluding wavecol) and transpose
    num_cols = len(data[0])
    outdata = []
    for col_idx in range(num_cols):
        if col_idx == wavecol:
            continue
        column = [row[col_idx] for row in data]
        outdata.append(column)

    return waves, outdata, out_headers

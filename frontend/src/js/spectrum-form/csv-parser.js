import { inferSchema, initParser } from "udsv"

/**
 * CSV/TSV Parser using uDSV
 *
 * Parses spectral data files with automatic delimiter and decimal format detection.
 * Supports European (semicolon/comma) and US (comma/period) formats.
 */

/**
 * Parse CSV/TSV text into structured data.
 *
 * @param {string} text - Raw file content
 * @returns {{ headers: string[], rows: number[][], delimiter: string, hasHeaders: boolean }}
 */
export function parseCSV(text) {
  if (!text || !text.trim()) {
    throw new Error("File is empty")
  }

  // Use udsv to infer schema (delimiter, headers, types)
  const schema = inferSchema(text)

  // Add custom parsers for European format (semicolon delimiter with comma decimals)
  if (schema.col === ";") {
    for (const col of schema.cols) {
      col.parse = (str) => {
        if (!str || typeof str !== "string") return str
        const trimmed = str.trim()
        // European: dots are thousands separators, commas are decimals
        const normalized = trimmed.replace(/\./g, "").replace(",", ".")
        const num = parseFloat(normalized)
        return Number.isNaN(num) ? str : num
      }
    }
  }

  const parser = initParser(schema)

  // Parse as typed objects (custom parsers only work with typed methods)
  const objRows = parser.typedObjs(text)

  if (!objRows || objRows.length === 0) {
    throw new Error("No data rows found in file")
  }

  // Extract headers from schema
  const headers = schema.cols.map((col) => col.name)

  // Validate: must have at least 2 columns
  if (headers.length < 2) {
    const detectedDelim =
      schema.col === "\t" ? "tab" : schema.col === "," ? "comma" : JSON.stringify(schema.col)
    throw new Error(
      `File must have at least 2 columns (wavelength and data). ` +
        `Only 1 column was detected (delimiter: ${detectedDelim}).`
    )
  }

  // Convert objects to arrays matching header order
  const rows = objRows.map((obj) => headers.map((header) => obj[header]))

  // Filter out rows with no valid numeric data
  const validRows = rows.filter((row) => row.some((v) => typeof v === "number" && !Number.isNaN(v)))

  if (validRows.length === 0) {
    throw new Error("No valid numeric data found in file")
  }

  // Validate: check that we have wavelength-like data in at least one column
  const hasWavelengthColumn = headers.some((_, colIndex) => {
    const colValues = validRows.map((row) => row[colIndex]).filter((v) => typeof v === "number")
    if (colValues.length < 3) return false
    // Check if values are in plausible wavelength range (150-1800nm)
    const inRange = colValues.filter((v) => v >= 150 && v <= 1800)
    return inRange.length >= colValues.length * 0.5 // At least 50% in range
  })

  if (!hasWavelengthColumn) {
    throw new Error(
      "Could not find a wavelength column. " +
        "Expected numeric values in the 150-1800nm range. " +
        "Check that your file format is correct."
    )
  }

  return {
    headers,
    rows: validRows,
    delimiter: schema.col,
    hasHeaders: true, // udsv treats first row as headers by default
  }
}

/**
 * Extract spectrum data from parsed CSV.
 *
 * @param {{ rows: number[][] }} parsed - Parsed CSV result
 * @param {number} waveColIndex - Index of wavelength column
 * @param {number} dataColIndex - Index of data column
 * @returns {number[][]} Array of [wavelength, value] pairs sorted by wavelength
 */
export function extractSpectrum(parsed, waveColIndex, dataColIndex) {
  const data = []

  for (const row of parsed.rows) {
    const wavelength = row[waveColIndex]
    const value = row[dataColIndex]

    if (Number.isNaN(wavelength) || Number.isNaN(value)) continue
    if (wavelength < 150 || wavelength > 1800) continue

    data.push([wavelength, value])
  }

  return data.sort((a, b) => a[0] - b[0])
}

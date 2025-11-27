/**
 * CSV/TSV Parser Module
 *
 * Parses text files with automatic delimiter and decimal format detection.
 * Supports both European (;/,) and US (,/.) formats.
 */

/**
 * Parse CSV/TSV text into a structured result.
 *
 * @param {string} text - Raw file content
 * @returns {{ headers: string[] | null, rows: number[][], delimiter: string, hasHeaders: boolean }}
 */
export function parseCSV(text) {
  // Normalize line endings
  text = text.replace(/\r\n/g, "\n").replace(/\r/g, "\n")

  // Auto-detect delimiter by counting occurrences in first few lines
  const sampleLines = text.split("\n").slice(0, 10).join("\n")
  const delimiter = detectDelimiter(sampleLines)
  const decimalChar = delimiter === ";" ? "," : "."

  // Split into lines and parse
  const lines = text.split("\n").filter((line) => line.trim())
  if (lines.length === 0) {
    throw new Error("File is empty")
  }

  // Parse all rows (as strings first)
  const rawRows = lines.map((line) => splitLine(line, delimiter))

  // Validate all rows have same column count
  const columnCounts = new Set(rawRows.map((r) => r.length))
  if (columnCounts.size > 1) {
    // Filter out rows with different lengths (likely trailing data)
    const mostCommon = [...columnCounts].reduce((a, b) =>
      rawRows.filter((r) => r.length === a).length > rawRows.filter((r) => r.length === b).length
        ? a
        : b
    )
    // Keep only rows with the most common length
    const filteredRows = rawRows.filter((r) => r.length === mostCommon)
    if (filteredRows.length < 2) {
      throw new Error("Inconsistent column counts in file")
    }
    rawRows.length = 0
    rawRows.push(...filteredRows)
  }

  // Detect if first row is headers
  const hasHeaders = detectHeaders(rawRows[0], decimalChar)

  // Extract headers and data rows
  let headers = null
  let dataRows = rawRows
  if (hasHeaders) {
    headers = rawRows[0]
    dataRows = rawRows.slice(1)
  } else {
    // Generate generic headers
    headers = rawRows[0].map((_, i) => `Column ${i + 1}`)
  }

  // Parse data rows to numbers
  const rows = dataRows.map((row) => row.map((cell) => parseNumber(cell, decimalChar)))

  // Filter out rows with all NaN
  const validRows = rows.filter((row) => row.some((v) => !Number.isNaN(v)))

  return {
    headers,
    rows: validRows,
    delimiter,
    hasHeaders,
  }
}

/**
 * Extract spectrum data from parsed CSV.
 *
 * @param {{ rows: number[][], headers: string[] }} parsed - Parsed CSV result
 * @param {number} waveColIndex - Index of wavelength column
 * @param {number} dataColIndex - Index of data column
 * @returns {number[][]} Array of [wavelength, value] pairs sorted by wavelength
 */
export function extractSpectrum(parsed, waveColIndex, dataColIndex) {
  const data = []

  for (const row of parsed.rows) {
    const wavelength = row[waveColIndex]
    const value = row[dataColIndex]

    // Skip if either value is NaN
    if (Number.isNaN(wavelength) || Number.isNaN(value)) {
      continue
    }

    // Skip if wavelength is outside reasonable range
    if (wavelength < 150 || wavelength > 1800) {
      continue
    }

    data.push([wavelength, value])
  }

  // Sort by wavelength
  data.sort((a, b) => a[0] - b[0])

  return data
}

/**
 * Detect the most likely delimiter in the text.
 */
function detectDelimiter(text) {
  const counts = {
    "\t": (text.match(/\t/g) || []).length,
    ";": (text.match(/;/g) || []).length,
    ",": (text.match(/,/g) || []).length,
  }

  // Tab has priority if present
  if (counts["\t"] > 0) {
    return "\t"
  }

  // Semicolon typically indicates European format
  if (counts[";"] > counts[","]) {
    return ";"
  }

  // Default to comma
  return ","
}

/**
 * Split a line by delimiter, handling quoted fields.
 */
function splitLine(line, delimiter) {
  const result = []
  let current = ""
  let inQuotes = false

  for (let i = 0; i < line.length; i++) {
    const char = line[i]

    if (char === '"') {
      inQuotes = !inQuotes
    } else if (char === delimiter && !inQuotes) {
      result.push(current.trim())
      current = ""
    } else {
      current += char
    }
  }
  result.push(current.trim())

  return result
}

/**
 * Detect if a row looks like headers (contains non-numeric values).
 */
function detectHeaders(row, decimalChar) {
  // If any cell can't be parsed as a number, it's likely a header
  for (const cell of row) {
    const trimmed = cell.trim()
    if (!trimmed) continue

    const parsed = parseNumber(trimmed, decimalChar)
    if (Number.isNaN(parsed)) {
      return true
    }
  }
  return false
}

/**
 * Parse a string to number, handling different decimal formats.
 */
function parseNumber(str, decimalChar) {
  if (typeof str !== "string") return NaN

  let normalized = str.trim()
  if (!normalized) return NaN

  // Handle European format (comma as decimal)
  if (decimalChar === ",") {
    // Remove dots (thousands separator)
    normalized = normalized.replace(/\./g, "")
    // Convert comma to dot (decimal)
    normalized = normalized.replace(",", ".")
  } else {
    // Remove commas (thousands separator in US format)
    normalized = normalized.replace(/,/g, "")
  }

  return parseFloat(normalized)
}

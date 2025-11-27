/**
 * CSV/TSV Parser
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
  const normalized = text.replace(/\r\n/g, "\n").replace(/\r/g, "\n")
  const lines = normalized.split("\n").filter((line) => line.trim())

  if (lines.length === 0) {
    throw new Error("File is empty")
  }

  const sampleText = lines.slice(0, 10).join("\n")
  const delimiter = detectDelimiter(sampleText)
  const decimalChar = delimiter === ";" ? "," : "."

  let rawRows = lines.map((line) => splitLine(line, delimiter))

  // Normalize to most common column count (handles trailing data)
  rawRows = normalizeRowLengths(rawRows)

  const hasHeaders = detectHeaders(rawRows[0], decimalChar)
  const headers = hasHeaders ? rawRows[0] : rawRows[0].map((_, i) => `Column ${i + 1}`)
  const dataRows = hasHeaders ? rawRows.slice(1) : rawRows

  const rows = dataRows
    .map((row) => row.map((cell) => parseNumber(cell, decimalChar)))
    .filter((row) => row.some((v) => !Number.isNaN(v)))

  return { headers, rows, delimiter, hasHeaders }
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

/**
 * Detect the most likely delimiter by counting occurrences.
 */
function detectDelimiter(text) {
  const counts = {
    "\t": (text.match(/\t/g) || []).length,
    ";": (text.match(/;/g) || []).length,
    ",": (text.match(/,/g) || []).length,
  }

  if (counts["\t"] > 0) return "\t"
  if (counts[";"] > counts[","]) return ";"
  return ","
}

/**
 * Split a line by delimiter, respecting quoted fields.
 */
function splitLine(line, delimiter) {
  const result = []
  let current = ""
  let inQuotes = false

  for (const char of line) {
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
 * Normalize rows to the most common column count.
 */
function normalizeRowLengths(rows) {
  const countFrequency = new Map()
  for (const row of rows) {
    countFrequency.set(row.length, (countFrequency.get(row.length) || 0) + 1)
  }

  let mostCommonLength = 0
  let maxFreq = 0
  for (const [length, freq] of countFrequency) {
    if (freq > maxFreq) {
      maxFreq = freq
      mostCommonLength = length
    }
  }

  const filtered = rows.filter((r) => r.length === mostCommonLength)
  if (filtered.length < 2) {
    throw new Error("Inconsistent column counts in file")
  }

  return filtered
}

/**
 * Detect if the first row contains headers (non-numeric values).
 */
function detectHeaders(row, decimalChar) {
  return row.some((cell) => {
    const trimmed = cell.trim()
    return trimmed && Number.isNaN(parseNumber(trimmed, decimalChar))
  })
}

/**
 * Parse a string to number, handling different decimal formats.
 */
function parseNumber(str, decimalChar) {
  if (typeof str !== "string") return NaN

  let normalized = str.trim()
  if (!normalized) return NaN

  if (decimalChar === ",") {
    // European: dots are thousands, commas are decimals
    normalized = normalized.replace(/\./g, "").replace(",", ".")
  } else {
    // US: commas are thousands
    normalized = normalized.replace(/,/g, "")
  }

  return parseFloat(normalized)
}

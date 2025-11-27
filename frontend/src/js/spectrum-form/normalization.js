/**
 * Normalization Module
 *
 * Client-side spectrum normalization and interpolation.
 * Ports the Python logic from proteins/util/spectra.py
 */

/**
 * Interpolate spectrum data to 1nm wavelength grid.
 *
 * @param {number[][]} data - Array of [wavelength, value] pairs
 * @returns {number[][]} Interpolated data with integer wavelengths at 1nm steps
 */
export function interpolateToOneNm(data) {
  if (data.length < 2) return data

  // Ensure data is sorted by wavelength
  const sorted = [...data].sort((a, b) => a[0] - b[0])

  // Get wavelength range (round to integers)
  const minWave = Math.ceil(sorted[0][0])
  const maxWave = Math.floor(sorted[sorted.length - 1][0])

  // Create output array with 1nm steps
  const result = []

  for (let wave = minWave; wave <= maxWave; wave++) {
    const y = linearInterpolate(sorted, wave)
    result.push([wave, y])
  }

  return result
}

/**
 * Linear interpolation at a given x value.
 */
function linearInterpolate(sortedData, x) {
  // Find surrounding points
  let i = 0
  while (i < sortedData.length - 1 && sortedData[i + 1][0] < x) {
    i++
  }

  // Edge cases
  if (i >= sortedData.length - 1) {
    return sortedData[sortedData.length - 1][1]
  }

  const [x0, y0] = sortedData[i]
  const [x1, y1] = sortedData[i + 1]

  // Avoid division by zero
  if (x1 === x0) return y0

  // Linear interpolation formula
  return y0 + ((y1 - y0) * (x - x0)) / (x1 - x0)
}

/**
 * Normalize spectrum to peak = 1.0
 *
 * @param {number[][]} data - Array of [wavelength, value] pairs (already interpolated)
 * @param {object} options - Normalization options
 * @param {number} [options.rangeMin] - Minimum wavelength for peak search
 * @param {number} [options.rangeMax] - Maximum wavelength for peak search
 * @returns {{ normalized: number[][], peakWave: number, peakValue: number }}
 */
export function normalizeSpectrum(data, options = {}) {
  if (data.length === 0) {
    return { normalized: [], peakWave: null, peakValue: 0 }
  }

  const { rangeMin, rangeMax } = options

  // Find peak within specified range (or entire spectrum)
  let peakIndex = 0
  let peakValue = -Infinity

  for (let i = 0; i < data.length; i++) {
    const [wave, value] = data[i]

    // Skip if outside range (check for both null and undefined)
    if (rangeMin != null && wave < rangeMin) continue
    if (rangeMax != null && wave > rangeMax) continue

    if (value > peakValue) {
      peakValue = value
      peakIndex = i
    }
  }

  const peakWave = data[peakIndex][0]

  // Normalize all values
  const normalized = data.map(([wave, value]) => {
    const normalizedValue = peakValue > 0 ? Math.max(value / peakValue, 0) : 0
    // Round to 4 decimal places (matches Python behavior)
    return [wave, Math.round(normalizedValue * 10000) / 10000]
  })

  return {
    normalized,
    peakWave,
    peakValue,
  }
}

/**
 * Normalize 2P spectrum using local maxima detection.
 *
 * For 2-photon spectra, we need to find the biologically relevant peak,
 * which may not be the absolute maximum (often the left tail is higher).
 *
 * @param {number[][]} data - Array of [wavelength, value] pairs
 * @param {object} options - Normalization options
 * @param {number} [options.rangeMin] - Minimum wavelength for peak search
 * @param {number} [options.rangeMax] - Maximum wavelength for peak search
 * @param {number} [options.order=100] - Number of points on each side for local max detection
 * @returns {{ normalized: number[][], peakWave: number, peakValue: number }}
 */
export function normalize2P(data, options = {}) {
  if (data.length === 0) {
    return { normalized: [], peakWave: null, peakValue: 0 }
  }

  const { rangeMin, rangeMax, order = 100 } = options

  // Extract y values
  const yValues = data.map(([_, y]) => y)

  // Find local maxima
  const localMaxIndices = findLocalMaxima(yValues, order)

  // Filter: skip first 10 points to avoid noise at spectrum start
  const validMaxima = localMaxIndices.filter((i) => i > 10)

  // Filter by range if specified (check for both null and undefined)
  const inRangeMaxima = validMaxima.filter((i) => {
    const wave = data[i][0]
    if (rangeMin != null && wave < rangeMin) return false
    if (rangeMax != null && wave > rangeMax) return false
    return true
  })

  // Find the highest local maximum
  let peakIndex = 0
  let peakValue = 0

  if (inRangeMaxima.length > 0) {
    for (const i of inRangeMaxima) {
      if (yValues[i] > peakValue) {
        peakValue = yValues[i]
        peakIndex = i
      }
    }
  } else {
    // Fallback: use absolute max in range
    let maxVal = -Infinity
    for (let i = 0; i < data.length; i++) {
      const [wave, value] = data[i]
      if (rangeMin != null && wave < rangeMin) continue
      if (rangeMax != null && wave > rangeMax) continue
      if (value > maxVal) {
        maxVal = value
        peakIndex = i
        peakValue = value
      }
    }
  }

  const peakWave = data[peakIndex][0]

  // Normalize all values
  const normalized = data.map(([wave, value]) => {
    const normalizedValue = peakValue > 0 ? Math.max(value / peakValue, 0) : 0
    return [wave, Math.round(normalizedValue * 10000) / 10000]
  })

  return {
    normalized,
    peakWave,
    peakValue,
  }
}

/**
 * Find local maxima in an array.
 *
 * A point is a local maximum if it's greater than all points
 * within 'order' distance on both sides.
 *
 * @param {number[]} values - Array of values
 * @param {number} order - Number of points on each side to compare
 * @returns {number[]} Indices of local maxima
 */
export function findLocalMaxima(values, order = 1) {
  const maxima = []

  for (let i = 0; i < values.length; i++) {
    let isMax = true

    // Check left neighbors
    const leftStart = Math.max(0, i - order)
    for (let j = leftStart; j < i; j++) {
      if (values[j] >= values[i]) {
        isMax = false
        break
      }
    }

    if (!isMax) continue

    // Check right neighbors
    const rightEnd = Math.min(values.length - 1, i + order)
    for (let j = i + 1; j <= rightEnd; j++) {
      if (values[j] >= values[i]) {
        isMax = false
        break
      }
    }

    if (isMax) {
      maxima.push(i)
    }
  }

  return maxima
}

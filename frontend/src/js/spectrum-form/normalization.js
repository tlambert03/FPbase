/**
 * Spectrum Normalization
 *
 * Client-side spectrum interpolation and normalization.
 * Ports logic from backend proteins/util/spectra.py
 */

/**
 * Interpolate spectrum data to 1nm wavelength grid.
 *
 * @param {number[][]} data - Array of [wavelength, value] pairs
 * @returns {number[][]} Interpolated data with integer wavelengths
 */
export function interpolateToOneNm(data) {
  if (data.length < 2) return data

  const sorted = [...data].sort((a, b) => a[0] - b[0])
  const minWave = Math.ceil(sorted[0][0])
  const maxWave = Math.floor(sorted[sorted.length - 1][0])
  const result = []

  for (let wave = minWave; wave <= maxWave; wave++) {
    result.push([wave, linearInterpolate(sorted, wave)])
  }

  return result
}

/**
 * Normalize spectrum to peak = 1.0
 *
 * @param {number[][]} data - Array of [wavelength, value] pairs
 * @param {object} options - Options
 * @param {number} [options.rangeMin] - Min wavelength for peak search
 * @param {number} [options.rangeMax] - Max wavelength for peak search
 * @returns {{ normalized: number[][], peakWave: number|null, peakValue: number }}
 */
export function normalizeSpectrum(data, options = {}) {
  if (data.length === 0) {
    return { normalized: [], peakWave: null, peakValue: 0 }
  }

  const { rangeMin, rangeMax } = options
  let peakIndex = 0
  let peakValue = -Infinity

  for (let i = 0; i < data.length; i++) {
    const [wave, value] = data[i]
    if (rangeMin != null && wave < rangeMin) continue
    if (rangeMax != null && wave > rangeMax) continue

    if (value > peakValue) {
      peakValue = value
      peakIndex = i
    }
  }

  const peakWave = data[peakIndex][0]
  const normalized = data.map(([wave, value]) => [
    wave,
    roundTo4(peakValue > 0 ? Math.max(value / peakValue, 0) : 0),
  ])

  return { normalized, peakWave, peakValue }
}

/**
 * Normalize 2P spectrum using local maxima detection.
 *
 * For 2-photon spectra, we find the biologically relevant peak which may
 * not be the absolute maximum (the left tail is often higher).
 *
 * @param {number[][]} data - Array of [wavelength, value] pairs
 * @param {object} options - Options
 * @param {number} [options.rangeMin] - Min wavelength for peak search
 * @param {number} [options.rangeMax] - Max wavelength for peak search
 * @param {number} [options.order=100] - Points on each side for local max detection
 * @returns {{ normalized: number[][], peakWave: number|null, peakValue: number }}
 */
export function normalize2P(data, options = {}) {
  if (data.length === 0) {
    return { normalized: [], peakWave: null, peakValue: 0 }
  }

  const { rangeMin, rangeMax, order = 100 } = options
  const yValues = data.map(([, y]) => y)
  const localMaxIndices = findLocalMaxima(yValues, order)

  // Filter: skip first 10 points (noise) and respect range
  const validMaxima = localMaxIndices.filter((i) => {
    if (i <= 10) return false
    const wave = data[i][0]
    if (rangeMin != null && wave < rangeMin) return false
    if (rangeMax != null && wave > rangeMax) return false
    return true
  })

  let peakIndex = 0
  let peakValue = 0

  if (validMaxima.length > 0) {
    for (const i of validMaxima) {
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
  const normalized = data.map(([wave, value]) => [
    wave,
    roundTo4(peakValue > 0 ? Math.max(value / peakValue, 0) : 0),
  ])

  return { normalized, peakWave, peakValue }
}

/**
 * Find local maxima in an array.
 *
 * A point is a local maximum if it's greater than all points
 * within 'order' distance on both sides.
 *
 * @param {number[]} values - Array of values
 * @param {number} order - Points on each side to compare
 * @returns {number[]} Indices of local maxima
 */
export function findLocalMaxima(values, order = 1) {
  const maxima = []

  for (let i = 0; i < values.length; i++) {
    let isMax = true

    // Check left neighbors
    for (let j = Math.max(0, i - order); j < i && isMax; j++) {
      if (values[j] >= values[i]) isMax = false
    }

    // Check right neighbors
    for (let j = i + 1; j <= Math.min(values.length - 1, i + order) && isMax; j++) {
      if (values[j] >= values[i]) isMax = false
    }

    if (isMax) maxima.push(i)
  }

  return maxima
}

// ============================================================================
// Helpers
// ============================================================================

function linearInterpolate(sortedData, x) {
  let i = 0
  while (i < sortedData.length - 1 && sortedData[i + 1][0] < x) {
    i++
  }

  if (i >= sortedData.length - 1) {
    return sortedData[sortedData.length - 1][1]
  }

  const [x0, y0] = sortedData[i]
  const [x1, y1] = sortedData[i + 1]

  if (x1 === x0) return y0
  return y0 + ((y1 - y0) * (x - x0)) / (x1 - x0)
}

function roundTo4(value) {
  return Math.round(value * 10000) / 10000
}

import type { ChartOptions } from "../types"

/**
 * Boolean chart options that can be serialized to URL
 * These are converted to 0/1 in URLs for backwards compatibility
 */
const BOOLEAN_CHART_OPTIONS = [
  "showY",
  "showX",
  "showGrid",
  "areaFill",
  "logScale",
  "scaleEC",
  "scaleQY",
  "shareTooltip",
] as const

/**
 * Partial state that can be serialized to/from URL parameters
 */
export interface SpectraURLState {
  activeSpectra?: string[]
  activeOverlaps?: string[]
  chartOptions?: Partial<ChartOptions>
  exNorm?: readonly [number, string] | null
}

/**
 * Parse URL search string into partial spectra state
 * Pure function - no side effects, no direct store access
 *
 * @param search - URL search string (e.g., "?s=17,18&showY=1")
 * @returns Parsed state object
 */
export function parseURLParams(search: string): SpectraURLState {
  const params = new URLSearchParams(search)

  // If no URL params, return empty object
  if (params.size === 0) {
    return {}
  }

  const result: SpectraURLState = {}

  // Parse spectra IDs (comma-separated string like "17,18")
  const spectraParam = params.get("s")
  if (spectraParam) {
    const spectraIds = spectraParam.split(",").filter(Boolean)
    if (spectraIds.length > 0) {
      result.activeSpectra = spectraIds
    }
  }

  // Parse overlap IDs (comma-separated string)
  const overlapParam = params.get("o")
  if (overlapParam) {
    const overlapIds = overlapParam.split(",").filter(Boolean)
    if (overlapIds.length > 0) {
      result.activeOverlaps = overlapIds
    }
  }

  // Parse excitation normalization (support both new 'ex' and legacy 'normWave'+'normID' formats)
  const exParam = params.get("ex")
  if (exParam) {
    const exIds = exParam.split(",").filter(Boolean)
    if (exIds.length === 2) {
      const wave = Number.parseFloat(exIds[0]!)
      const id = exIds[1]!
      if (!Number.isNaN(wave) && id) {
        result.exNorm = [wave, id] as const
      }
    }
  } else {
    // Legacy format: normWave=452&normID=%24cl0
    const normWave = params.get("normWave")
    const normID = params.get("normID")
    if (normWave && normID) {
      const wave = Number.parseFloat(normWave)
      if (!Number.isNaN(wave)) {
        result.exNorm = [wave, normID] as const
      }
    }
  }

  // Parse chart options
  const chartOptions: Partial<ChartOptions> = {}
  let hasChartOptions = false

  // Boolean options
  for (const option of BOOLEAN_CHART_OPTIONS) {
    const value = params.get(option)
    if (value !== null) {
      // URL params are strings, convert to boolean
      chartOptions[option] = value === "true" || value === "1"
      hasChartOptions = true
    }
  }

  // String options
  const palette = params.get("palette")
  if (palette) {
    chartOptions.palette = palette
    hasChartOptions = true
  }

  // Extremes (support both 'min'+'max' and legacy 'xMin'+'xMax' formats)
  const minParam = params.get("min") || params.get("xMin")
  const maxParam = params.get("max") || params.get("xMax")
  if (minParam !== null && maxParam !== null) {
    const min = Number.parseFloat(minParam)
    const max = Number.parseFloat(maxParam)
    if (!Number.isNaN(min) && !Number.isNaN(max)) {
      chartOptions.extremes = [min, max]
      hasChartOptions = true
    }
  }

  if (hasChartOptions) {
    result.chartOptions = chartOptions
  }

  return result
}

/**
 * Serialize partial spectra state into URL search string
 * Pure function - no side effects, no direct store access
 * Maintains backwards compatibility with 6+ year old URL format
 *
 * @param state - Partial state to serialize
 * @returns URL search string (without leading '?')
 */
export function serializeURLParams(state: SpectraURLState): string {
  const params = new URLSearchParams()

  // Add spectra IDs (comma-separated)
  if (state.activeSpectra && state.activeSpectra.length > 0) {
    params.set("s", state.activeSpectra.join(","))
  }

  // Add overlap IDs (comma-separated)
  if (state.activeOverlaps && state.activeOverlaps.length > 0) {
    params.set("o", state.activeOverlaps.join(","))
  }

  // Add chart options
  if (state.chartOptions) {
    // Boolean options - convert to 0/1 for backwards compatibility
    for (const key of BOOLEAN_CHART_OPTIONS) {
      if (state.chartOptions[key] !== undefined) {
        params.set(key, state.chartOptions[key] ? "1" : "0")
      }
    }

    // String options
    if (state.chartOptions.palette) {
      params.set("palette", state.chartOptions.palette)
    }

    // Extremes - use xMin/xMax for backwards compatibility
    if (state.chartOptions.extremes && Array.isArray(state.chartOptions.extremes)) {
      const [xMin, xMax] = state.chartOptions.extremes
      if (xMin !== null && xMin !== undefined) {
        params.set("xMin", String(xMin))
      }
      if (xMax !== null && xMax !== undefined) {
        params.set("xMax", String(xMax))
      }
    }
  }

  // Add excitation normalization - use normWave/normID for backwards compatibility
  if (state.exNorm && state.exNorm[0] !== null && state.exNorm[1] !== null) {
    params.set("normWave", String(state.exNorm[0]))
    params.set("normID", String(state.exNorm[1]))
  }

  return params.toString()
}

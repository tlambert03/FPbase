import { defaultChartOptions } from "../defaults"
import type { ChartOptions, CustomFilterParams, CustomLaserParams } from "../types"

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
  hiddenSpectra?: string[]
  chartOptions?: Partial<ChartOptions>
  exNorm?: readonly [number, string] | null
  customFilters?: Record<string, CustomFilterParams>
  customLasers?: Record<string, CustomLaserParams>
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

  // Parse spectra IDs (comma-separated string like "17,18,$cf0_BP_450_50_90,$cl1_488")
  // For custom spectra, extract stable IDs and parameters
  const spectraParam = params.get("s")
  if (spectraParam) {
    const spectraIds: string[] = []
    const customFilters: Record<string, CustomFilterParams> = {}
    const customLasers: Record<string, CustomLaserParams> = {}

    for (const id of spectraParam.split(",").filter(Boolean)) {
      if (id.startsWith("$cf")) {
        // Parse custom filter: $cf0_BP_450_50_90
        const parts = id.split("_")
        const stableId = parts[0]
        if (stableId && parts.length >= 5 && parts[1] && parts[2] && parts[3] && parts[4]) {
          const type = parts[1].toUpperCase() as "BP" | "LP" | "SP"
          const center = Number.parseInt(parts[2], 10)
          const width = Number.parseInt(parts[3], 10)
          const transmission = Number.parseInt(parts[4], 10)

          spectraIds.push(stableId)
          customFilters[stableId] = { type, center, width, transmission }
        } else {
          // Malformed, just add the ID
          spectraIds.push(id)
        }
      } else if (id.startsWith("$cl")) {
        // Parse custom laser: $cl1_488
        const parts = id.split("_")
        const stableId = parts[0]
        if (stableId && parts.length >= 2 && parts[1]) {
          const wavelength = Number.parseInt(parts[1], 10)

          spectraIds.push(stableId)
          customLasers[stableId] = { wavelength }
        } else {
          // Malformed, just add the ID
          spectraIds.push(id)
        }
      } else {
        // Regular spectrum ID
        spectraIds.push(id)
      }
    }

    if (spectraIds.length > 0) {
      result.activeSpectra = spectraIds
    }
    if (Object.keys(customFilters).length > 0) {
      result.customFilters = customFilters
    }
    if (Object.keys(customLasers).length > 0) {
      result.customLasers = customLasers
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

  // Parse hidden spectra IDs (comma-separated string)
  const hiddenParam = params.get("h")
  if (hiddenParam) {
    const hiddenIds = hiddenParam.split(",").filter(Boolean)
    if (hiddenIds.length > 0) {
      result.hiddenSpectra = hiddenIds
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
  // Allow individual extremes to be set (e.g., only min or only max)
  const minParam = params.get("min") || params.get("xMin")
  const maxParam = params.get("max") || params.get("xMax")
  if (minParam !== null || maxParam !== null) {
    const min = minParam !== null ? Number.parseFloat(minParam) : null
    const max = maxParam !== null ? Number.parseFloat(maxParam) : null

    // Only set extremes if at least one valid number was parsed
    if ((min !== null && !Number.isNaN(min)) || (max !== null && !Number.isNaN(max))) {
      chartOptions.extremes = [
        min !== null && !Number.isNaN(min) ? min : null,
        max !== null && !Number.isNaN(max) ? max : null,
      ]
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
  // For custom spectra, serialize with parameters for backward compatibility
  if (state.activeSpectra && state.activeSpectra.length > 0) {
    const serializedIds = state.activeSpectra.map((id) => {
      if (id.startsWith("$cf") && state.customFilters?.[id]) {
        // Serialize custom filter: $cf0 + params -> $cf0_BP_450_50_90
        const filter = state.customFilters[id]
        return `${id}_${filter.type}_${filter.center}_${filter.width}_${filter.transmission}`
      }
      if (id.startsWith("$cl") && state.customLasers?.[id]) {
        // Serialize custom laser: $cl1 + params -> $cl1_488
        const laser = state.customLasers[id]
        return `${id}_${laser.wavelength}`
      }
      // Regular spectrum ID
      return id
    })
    params.set("s", serializedIds.join(","))
  }

  // Add overlap IDs (comma-separated)
  if (state.activeOverlaps && state.activeOverlaps.length > 0) {
    params.set("o", state.activeOverlaps.join(","))
  }

  // Add hidden spectra IDs (comma-separated)
  if (state.hiddenSpectra && state.hiddenSpectra.length > 0) {
    params.set("h", state.hiddenSpectra.join(","))
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

/**
 * Canonicalize a partial state by filling in defaults for missing fields
 * Converts various partial representations into one standard complete form
 *
 * @param state - Partial state to canonicalize
 * @returns Canonical state with all fields present
 */
export function canonicalizeURLState(state: SpectraURLState): SpectraURLState {
  return {
    activeSpectra: state.activeSpectra ?? [],
    activeOverlaps: state.activeOverlaps ?? [],
    hiddenSpectra: state.hiddenSpectra ?? [],
    chartOptions: state.chartOptions
      ? { ...defaultChartOptions, ...state.chartOptions }
      : defaultChartOptions,
    exNorm: state.exNorm ?? null,
    customFilters: state.customFilters ?? {},
    customLasers: state.customLasers ?? {},
  }
}

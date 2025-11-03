import qs from "qs"
import type { ChartOptions, SpectraStore } from "../types"

interface URLParams extends Partial<ChartOptions> {
  s?: string | string[] // Spectra IDs
  o?: string | string[] // Overlap IDs
  ex?: string | string[] // Excitation normalization (new format)
  normWave?: string // Excitation normalization wavelength (legacy)
  normID?: string // Excitation normalization ID (legacy)
  min?: string // Range minimum
  max?: string // Range maximum
  xMin?: string // Range minimum (legacy)
  xMax?: string // Range maximum (legacy)
  [key: string]: string | string[] | boolean | number | [number, number] | null | undefined
}

// Parse URL parameters and update store
export function syncURLToStore(store: SpectraStore) {
  const params = qs.parse(window.location.search, {
    ignoreQueryPrefix: true,
  }) as URLParams

  // If no URL params, don't override store (use sessionStorage)
  if (Object.keys(params).length === 0) {
    return
  }

  // Parse spectra IDs (comma-separated string like "17,18" or array)
  if (params.s) {
    let spectraIds: string[]
    if (Array.isArray(params.s)) {
      spectraIds = params.s
    } else {
      // Split comma-separated string into array
      spectraIds = params.s.split(",")
    }
    store.setActiveSpectra(spectraIds.filter(Boolean))
  }

  // Parse overlap IDs (comma-separated string or array)
  if (params.o) {
    let overlapIds: string[]
    if (Array.isArray(params.o)) {
      overlapIds = params.o
    } else {
      // Split comma-separated string into array
      overlapIds = params.o.split(",")
    }
    store.setActiveOverlaps(overlapIds.filter(Boolean))
  }

  // Parse excitation normalization (support both new 'ex' and legacy 'normWave'+'normID' formats)
  if (params.ex) {
    let exIds: string[]
    if (Array.isArray(params.ex)) {
      exIds = params.ex
    } else {
      // Split comma-separated string into array
      exIds = params.ex.split(",")
    }
    const filtered = exIds.filter(Boolean)
    if (filtered.length === 2) {
      const wave = Number.parseFloat(filtered[0]!)
      const id = filtered[1]!
      if (!Number.isNaN(wave) && id) {
        store.setExNorm([wave, id] as const)
      }
    }
  } else if (params.normWave && params.normID) {
    // Legacy format: normWave=452&normID=%24cl0
    const wave = Number.parseFloat(params.normWave as string)
    const id = params.normID as string
    if (!Number.isNaN(wave) && id) {
      store.setExNorm([wave, id] as const)
    }
  }

  // Parse chart options
  const chartOptions: Partial<ChartOptions> = {}
  let hasChartOptions = false

  // Boolean options
  const booleanOptions = [
    "showY",
    "showX",
    "showGrid",
    "areaFill",
    "logScale",
    "scaleEC",
    "scaleQY",
    "shareTooltip",
  ] as const
  for (const option of booleanOptions) {
    const value = params[option]
    if (value !== undefined) {
      // URL params come as strings, convert to boolean
      chartOptions[option] = typeof value === "boolean" ? value : value === "true" || value === "1"
      hasChartOptions = true
    }
  }

  // String options
  if (params.palette) {
    chartOptions.palette = params.palette as string
    hasChartOptions = true
  }

  // Extremes (support both 'min'+'max' and legacy 'xMin'+'xMax' formats)
  const minParam = params.min || params.xMin
  const maxParam = params.max || params.xMax
  if (minParam !== undefined && maxParam !== undefined) {
    const min = Number.parseFloat(minParam as string)
    const max = Number.parseFloat(maxParam as string)
    if (!Number.isNaN(min) && !Number.isNaN(max)) {
      chartOptions.extremes = [min, max]
      hasChartOptions = true
    }
  }

  if (hasChartOptions) {
    store.updateChartOptions(chartOptions)
  }

  // Mark that store was initialized from URL
  store.setUrlInitialized(true)
}

// Hook to initialize URL sync on mount
export function useURLSync() {
  // This will be called once on component mount
  // We don't need to import useSpectraStore here since it will be
  // called from components that already have access to the store
  return { syncURLToStore }
}

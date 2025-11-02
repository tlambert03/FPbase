import qs from "qs"
import type { ChartOptions, SpectraStore } from "../types"

interface URLParams extends Partial<ChartOptions> {
  s?: string | string[] // Spectra IDs
  o?: string | string[] // Overlap IDs
  ex?: string | string[] // Excitation normalization
  min?: string // Range minimum
  max?: string // Range maximum
  [key: string]: string | string[] | boolean | number | null | undefined
}

// Parse URL parameters and update store
export function syncURLToStore(store: SpectraStore) {
  const params = qs.parse(window.location.search, {
    ignoreQueryPrefix: true,
  }) as URLParams

  // Parse spectra IDs
  if (params.s) {
    const spectraIds = Array.isArray(params.s) ? params.s : [params.s]
    store.setActiveSpectra(spectraIds.filter(Boolean))
  }

  // Parse overlap IDs
  if (params.o) {
    const overlapIds = Array.isArray(params.o) ? params.o : [params.o]
    store.setActiveOverlaps(overlapIds.filter(Boolean))
  }

  // Parse excitation normalization
  if (params.ex) {
    const exIds = Array.isArray(params.ex) ? params.ex : [params.ex]
    store.setExNorm(exIds.filter(Boolean))
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
  ]
  for (const option of booleanOptions) {
    if (params[option] !== undefined) {
      chartOptions[option] = params[option] === "true" || params[option] === "1"
      hasChartOptions = true
    }
  }

  // String options
  if (params.palette) {
    chartOptions.palette = params.palette
    hasChartOptions = true
  }

  // Extremes (min/max range)
  if (params.min !== undefined && params.max !== undefined) {
    const min = Number.parseFloat(params.min as string)
    const max = Number.parseFloat(params.max as string)
    if (!Number.isNaN(min) && !Number.isNaN(max)) {
      chartOptions.extremes = [min, max]
      hasChartOptions = true
    }
  }

  if (hasChartOptions) {
    store.updateChartOptions(chartOptions)
  }
}

// Hook to initialize URL sync on mount
export function useURLSync() {
  // This will be called once on component mount
  // We don't need to import useSpectraStore here since it will be
  // called from components that already have access to the store
  return { syncURLToStore }
}

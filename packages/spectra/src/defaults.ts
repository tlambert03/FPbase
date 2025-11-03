/**
 * Default chart options and initial state
 */

import type { ChartOptions, ExNorm, SpectrumSubtype } from "./types"

export const defaultChartOptions: ChartOptions = {
  showY: false,
  showX: true,
  showGrid: false,
  logScale: false,
  scaleEC: false,
  scaleQY: false,
  shareTooltip: true,
  areaFill: true,
  palette: "wavelength",
  extremes: null,
}

export interface Defaults {
  activeSpectra: string[]
  activeOverlaps: string[]
  chartOptions: ChartOptions
  exNorm: ExNorm
  excludeSubtypes: SpectrumSubtype[]
}

export const defaults: Defaults = {
  activeSpectra: [],
  activeOverlaps: [],
  chartOptions: defaultChartOptions,
  exNorm: null,
  excludeSubtypes: ["2P"], // Exclude two-photon spectra by default
}

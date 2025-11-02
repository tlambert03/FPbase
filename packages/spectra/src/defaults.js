/**
 * Default chart options and initial state
 */

export const defaultChartOptions = {
  showY: false,
  showX: true,
  showGrid: false,
  logScale: false,
  scaleEC: false,
  scaleQY: false,
  shareTooltip: true,
  areaFill: true,
  palette: "wavelength",
  extremes: [null, null],
}

export const defaults = {
  activeSpectra: [],
  activeOverlaps: [],
  chartOptions: defaultChartOptions,
  exNorm: null, // [wavelength, spectrumId] tuple
  excludeSubtypes: ["2P"],
}

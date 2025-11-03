/**
 * Generate shareable URL from current spectra viewer state
 * Maintains backwards compatibility with 6+ year old URL format
 */
function stateToUrl(activeSpectra, chartOptions, exNorm) {
  const params = new URLSearchParams()

  // Add spectra IDs (comma-separated)
  if (activeSpectra && activeSpectra.length > 0) {
    params.set("s", activeSpectra.join(","))
  }

  // Add chart options
  if (chartOptions) {
    // Boolean options - convert to 0/1 for backwards compatibility
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
    for (const key of booleanOptions) {
      if (chartOptions[key] !== undefined) {
        params.set(key, chartOptions[key] ? "1" : "0")
      }
    }

    // String options
    if (chartOptions.palette) {
      params.set("palette", chartOptions.palette)
    }

    // Extremes - use xMin/xMax for backwards compatibility
    if (chartOptions.extremes && Array.isArray(chartOptions.extremes)) {
      const [xMin, xMax] = chartOptions.extremes
      if (xMin !== null && xMin !== undefined) {
        params.set("xMin", String(xMin))
      }
      if (xMax !== null && xMax !== undefined) {
        params.set("xMax", String(xMax))
      }
    }
  }

  // Add excitation normalization - use normWave/normID for backwards compatibility
  if (exNorm && exNorm[0] !== null && exNorm[1] !== null) {
    params.set("normWave", String(exNorm[0]))
    params.set("normID", String(exNorm[1]))
  }

  // Build full URL
  const { origin, pathname } = window.location
  const queryString = params.toString()
  return queryString ? `${origin}${pathname}?${queryString}` : ""
}

export default stateToUrl

/**
 * Spectrum Chart
 *
 * Highcharts-based spectrum preview with click-to-set-peak support.
 */

import Highcharts from "highcharts"

/**
 * Create a spectrum preview chart.
 *
 * @param {string|HTMLElement} container - Container element or ID
 * @param {number[][]} data - Array of [wavelength, value] pairs
 * @param {object} options - Chart options
 * @param {string} [options.name="Spectrum"] - Series name
 * @param {string} [options.color="#0d6efd"] - Series color
 * @param {function} [options.onClick] - Callback when chart is clicked: (wavelength) => void
 * @returns {object} Chart controller with update methods
 */
export function createSpectrumChart(container, data, options = {}) {
  const { name = "Spectrum", color = "#0d6efd", onClick, normalized = true } = options

  const containerId = typeof container === "string" ? container : container.id
  const containerEl = typeof container === "string" ? document.getElementById(container) : container

  if (!containerEl) {
    throw new Error(`Container not found: ${container}`)
  }

  const wavelengths = data.map(([x]) => x)
  const minWave = Math.min(...wavelengths)
  const maxWave = Math.max(...wavelengths)

  const handleClick = (e) => {
    if (onClick) {
      const wavelength = e.xAxis?.[0]?.value ?? e.point?.x
      if (wavelength != null) onClick(wavelength)
    }
  }

  // Calculate y-axis config based on normalized flag
  const yAxisConfig = getYAxisConfig(data, normalized)

  const chart = Highcharts.chart(containerId, {
    chart: {
      type: "area",
      height: 300,
      backgroundColor: "#ffffff",
      animation: { duration: 200 },
      spacing: [10, 10, 25, 10],
      events: { click: handleClick },
    },
    title: { text: null },
    xAxis: {
      title: { text: "Wavelength (nm)" },
      min: minWave - 10,
      max: maxWave + 10,
      tickLength: 3,
      gridLineWidth: 0,
      plotLines: [],
    },
    yAxis: yAxisConfig,
    tooltip: {
      shared: true,
      headerFormat: "<b>{point.x} nm</b><br/>",
      pointFormat: "{series.name}: {point.y}",
    },
    legend: { enabled: false },
    plotOptions: {
      area: {
        animation: false,
        marker: { enabled: false },
        lineWidth: 2,
        fillOpacity: 0.3,
        threshold: 0,
        events: { click: (e) => onClick?.(e.point.x) },
      },
    },
    series: [{ name, data: data.map(([x, y]) => [x, y]), color }],
    credits: { enabled: false },
    accessibility: { enabled: false },
  })

  return {
    chart,

    updateData(newData, seriesName = name) {
      chart.series[0].setData(
        newData.map(([x, y]) => [x, y]),
        true
      )
      chart.series[0].update({ name: seriesName }, false)

      const newWavelengths = newData.map(([x]) => x)
      chart.xAxis[0].setExtremes(Math.min(...newWavelengths) - 10, Math.max(...newWavelengths) + 10)
    },

    updateYAxis(newData, isNormalized) {
      const yAxisConfig = getYAxisConfig(newData, isNormalized)
      chart.yAxis[0].update(yAxisConfig, true)
    },

    setPeakMarker(wavelength) {
      chart.xAxis[0].removePlotLine("peak-marker")

      if (wavelength !== null) {
        chart.xAxis[0].addPlotLine({
          id: "peak-marker",
          value: wavelength,
          color: "#36a83b",
          width: 2,
          dashStyle: "Dash",
          label: {
            text: `Peak: ${wavelength} nm`,
            style: { color: "#36a83b", fontWeight: "bold" },
            rotation: 0,
            y: 15,
          },
          zIndex: 5,
        })
      }
    },

    clearAnnotations() {
      chart.xAxis[0].removePlotLine("peak-marker")
    },

    destroy() {
      chart?.destroy()
    },
  }
}

/**
 * Get y-axis configuration based on data and normalization state.
 *
 * @param {number[][]} data - Array of [wavelength, value] pairs
 * @param {boolean} normalized - Whether data is normalized
 * @returns {object} Highcharts y-axis configuration
 */
function getYAxisConfig(data, normalized) {
  if (normalized) {
    return {
      title: { text: "Relative Intensity" },
      min: 0,
      max: 1.05,
      labels: {
        formatter() {
          return `${Math.round(this.value * 100)}%`
        },
      },
      gridLineWidth: 1,
      gridLineColor: "#e0e0e0",
    }
  }

  // For absolute values, auto-scale based on data
  const yValues = data.map(([, y]) => y)
  const maxY = Math.max(...yValues)
  const minY = Math.min(...yValues)
  const range = maxY - minY
  const padding = range * 0.05 // 5% padding

  return {
    title: { text: "Intensity" },
    min: Math.max(0, minY - padding),
    max: maxY + padding,
    labels: {
      formatter() {
        return this.value.toFixed(2)
      },
    },
    gridLineWidth: 1,
    gridLineColor: "#e0e0e0",
  }
}

const FONTS =
  'Roboto, -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";'

const toolTipRow = entry => {
  return (
    `<tr><td><span style="color:${entry.series.color}">` +
    `&#9673; ` +
    `</span>${entry.series.name}${
      entry.series.userOptions.scaleEC
        ? ' <span style="font-size: 0.7rem; font-style: italic">(EC)</span>'
        : ""
    }${
      entry.series.userOptions.scaleQY
        ? ' <span style="font-size: 0.7rem; font-style: italic">(QY)</span>'
        : ""
    }:</td><td style="text-align: right; font-weight: bold"> ${(entry.series
      .userOptions.scaleEC
      ? Math.round(entry.y)
      : (Math.round(100 * entry.y) / 100).toFixed(2)
    ).toLocaleString()}</td></tr>`
  )
}
const DEFAULT_OPTIONS = {
  plotOptions: {
    series: {
      boostThreshold: 800,
      events: {
        mouseOver: function({ target: { xAxis } }, b) {
          const el = document.getElementById("zoom-info")
          if (el) el.style.display = "block"
        },
        mouseOut: function() {
          const el = document.getElementById("zoom-info")
          if (el) el.style.display = "none"
        },
      },
      animation: false,
      lineWidth: 0.5,
      marker: {
        enabled: false,
        symbol: "circle",
      },
    },
  },
  boost: {
    useGPUTranslations: true,
    seriesThreshold: 4,
    enabled: true,
    debug: {
      showSkipSummary: true,
      timeSeriesProcessing: true,
      timeRendering: true,
    },
  },
  chart: {
    type: "areaspline",
    zoomType: "x",
    panning: true,
    panKey: "shift",
    animation: { duration: 50 },
    resetZoomButton: {
      theme: {
        fill: "RGBA(90, 90, 90, 0.35)",
        style: { color: "#FFFFFF", fontSize: "0.7rem", fontWeight: "bold" },
        r: 3,
        stroke: "RGBA(90, 90, 90, 0.1)",
        states: {
          hover: {
            fill: "RGBA(18, 75, 51, 0.8)",
          },
        },
      },
    },
  },
  yAxis: {
    gridLineWidth: 1,
    title: false,
    labels: {
      enabled: false,
      formatter: function() {
        if (this.value !== 0) {
          return this.axis.defaultLabelFormatter.call(this)
        }
        return null
      },
    },
    maxPadding: 0.01,
  },
  xAxis: {
    gridLineWidth: 1,
    tickLength: 0,
    labels: {
      y: 15,
      enabled: true,
    },
    crosshair: true,
    events: {
      afterSetExtremes: function({ userMin, userMax }) {
        const el = document.getElementById("zoom-info")
        if (el) {
          if (userMin || userMax) {
            if (window.USER_IS_TOUCHING) {
              el.innerHTML = "pinch to zoom, two-finger drag to pan"
            } else {
              el.innerHTML =
                "click and drag to zoom, shift-click and drag to pan"
            }
          } else if (window.USER_IS_TOUCHING) {
            el.innerHTML = "pinch to zoom"
          } else {
            el.innerHTML = "click and drag to zoom"
          }
        }
      },
    },
  },
  navigation: {
    buttonOptions: {
      theme: {
        "stroke-width": 0,
        opacity: 0.35,
        states: {
          hover: {
            fill: null,
            opacity: 0.9,
          },
          select: {
            fill: "#EEEEEE",
            opacity: 0.9,
          },
        },
      },
    },
    menuItemStyle: {
      color: "#333",
      fontFamily: FONTS,
      fontSize: "0.7rem",
    },
  },
  exporting: {
    filename: "FPbase_Spectra.csv",
    sourceWidth: 1200,
    scale: 1,
    csv: {},
    chartOptions: {
      chart: {
        height: this && this.chartHeight,
      },
      title: false,
    },

    buttons: {
      contextButton: {
        enabled: false,
        menuItems: [
          "downloadPNG",
          "downloadPDF",
          "downloadSVG",
          "separator",
          "downloadCSV",
          "printChart",
          // "separator",
          // "reset"
          // "openInCloud"
          // "viewData"
        ],
      },
    },
  },
  legend: {
    verticalAlign: "top",
    align: "right",
    labelFormatter: function() {
      let { name } = this
      if (this.chart.chartWidth < 800) {
        name = name.replace("Chroma", "Chr").replace("Semrock", "Sem")
      }
      if (+this.chart.chartWidth < 500) {
        name = name.length > 18 ? `${name.slice(0, 18)}...` : name
      }
      return name
    },
    itemStyle: {
      fontWeight: 600,
      fontSize: "11px",
      fontFamily: FONTS,
    },
  },
  tooltip: {
    useHTML: true,
    backgroundColor: "#FFFFFFCC",
    borderColor: "#999",
    crosshairs: true,
    shared: true,
    hideDelay: 150,
    valueDecimals: 3,
    formatter: function(tooltip) {
      let tooltipHtml = "<table class='spectrum-tooltip'>"
      tooltipHtml += `${"<tr><td></td>" +
        "<td style='text-align: right; line-height: 1.1rem; font-size: 0.75rem; border-bottom: 1px solid #ccc;'>"}${
        this.x
      }nm</td></tr>`

      if (this.point) {
        tooltipHtml += toolTipRow(this.point)
      } else {
        this.points.forEach(function(entry) {
          tooltipHtml += toolTipRow(entry)
        })
      }

      tooltipHtml += "</table>"
      return tooltipHtml
    },
    positioner(labelWidth, labelHeight, point) {
      const chartwidth = this.chart.chartWidth
      const yAx2 = this.chart.get("yAx2")
      const rightPad = (yAx2 && yAx2.axisTitleMargin) || 0
      const y = Math.min(
        Math.max(point.plotY, 50),
        this.chart.chartHeight - labelHeight - 40
      )
      const t = 10 + point.plotX + labelWidth / 3
      const x =
        t + labelWidth < chartwidth - rightPad
          ? t
          : point.plotX - labelWidth - 20
      return { x, y }
    },
    shadow: false,
    style: {
      color: "#333",
      fontFamily: FONTS,
      fontSize: "0.8rem",
    },
  },
}

export default DEFAULT_OPTIONS

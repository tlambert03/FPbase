const FONTS =
  'Roboto, -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";'

const toolTipRow = entry => {
  return (
    '<tr><td><span style="color:' +
    entry.series.color +
    '">' +
    "&#9673; " +
    "</span>" +
    entry.series.name +
    (entry.series.userOptions.scaleEC
      ? ' <span style="font-size: 0.7rem; font-style: italic">(EC)</span>'
      : "") +
    (entry.series.userOptions.scaleQY
      ? ' <span style="font-size: 0.7rem; font-style: italic">(QY)</span>'
      : "") +
    ':</td><td style="text-align: right; font-weight: bold"> ' +
    (entry.series.userOptions.scaleEC
      ? Math.round(entry.y)
      : (Math.round(100 * entry.y) / 100).toFixed(2)
    ).toLocaleString() +
    "</td></tr>"
  )
}
const DEFAULT_OPTIONS = {
  plotOptions: {
    series: {
      animation: false,
      lineWidth: 0.5,
      marker: {
        enabled: false,
        symbol: "circle"
      }
    }
  },
  chart: {
    type: "areaspline",
    zoomType: "x",
    panning: true,
    panKey: "shift",
    animation: { duration: 50 }
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
      }
    },
    maxPadding: 0.01
  },
  xAxis: {
    gridLineWidth: 1,
    tickLength: 0,
    labels: {
      y: 15,
      enabled: true
    },
    crosshair: true
  },
  navigation: {
    buttonOptions: {
      theme: {
        "stroke-width": 0,
        opacity: 0.35,
        states: {
          hover: {
            fill: null,
            opacity: 0.9
          },
          select: {
            fill: "#EEEEEE",
            opacity: 0.9
          }
        }
      }
    },
    menuItemStyle: {
      color: "#333",
      fontFamily: FONTS,
      fontSize: "0.7rem"
    }
  },
  exporting: {
    filename: 'FPbase_Spectra.csv',
    sourceWidth: 1200,
    scale: 1,
    csv: {},
    chartOptions: {
      chart: {
        height: this && this.chartHeight
      },
      title: false
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
          "printChart"
          //"separator",
          //"reset"
          // "openInCloud"
          // "viewData"
        ]
      }
    }
  },
  legend: {
    verticalAlign: "top",
//    align: "left",
    itemStyle: {
      fontWeight: 600,
      fontSize: "11px",
      fontFamily: FONTS
    }
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
      let tooltip_html = "<table class='spectrum-tooltip'>"
      tooltip_html +=
        "<tr><td></td>" +
        "<td style='text-align: right; line-height: 1.1rem; font-size: 0.75rem; border-bottom: 1px solid #ccc;'>" +
        this.x +
        "nm</td></tr>"

      if (this.point) {
        tooltip_html += toolTipRow(this.point)
      } else {
        this.points.forEach(function(entry) {
          tooltip_html += toolTipRow(entry)
        })
      }

      tooltip_html += "</table>"
      return tooltip_html
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
      fontSize: "0.8rem"
    }
  }
}

export default DEFAULT_OPTIONS

import React, { useRef } from "react";
import Highcharts from "highcharts";
import HighchartsReact from "highcharts-react-official";

require("highcharts/modules/exporting")(Highcharts);

const SPECTRUM_CHARTOPTS = {
  chart: {
    zoomType: "x",
    type: "areaspline",
    animation: { duration: 200 },
    panKey: "meta",
    panning: true,
    events: {
      click: e => console.log(e)
    }
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
    }
  },
  credits: false,
  legend: {
    verticalAlign: "top",
    align: "right",
    x: -20
  },
  tooltip: {
    useHTML: true,
    backgroundColor: "#FFFFFFCC",
    borderColor: "#999",
    crosshairs: true,
    shared: true,
    valueDecimals: 3,
    positioner(labelWidth, labelHeight, point) {
      const chartwidth = this.chart.chartWidth;
      const y = Math.min(
        Math.max(point.plotY, 60),
        this.chart.chartHeight - labelHeight - 40
      );
      if (40 + point.plotX + labelWidth < chartwidth) {
        return {
          x: point.plotX + 35,
          y
        };
      }
      return {
        x: point.plotX - labelWidth - 20,
        y
      };
    },
    shadow: false,
    style: {
      color: "#333",
      fontFamily:
        '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";',
      fontSize: "0.8rem"
    }
  },
  series: [],
  title: false,
  xAxis: {
    tickLength: 0,
    labels: {
      y: 15
    }
  },
  yAxis: {
    min: 0,
    max: 1,
    gridLineWidth: 0,
    labels: false,
    title: false
  },
  plotOptions: {
    series: {
      fillOpacity: 0.45,
      animation: false,
      lineWidth: 0,
      marker: {
        enabled: false,
        symbol: "circle"
      },
      states: {
        hover: {
          halo: false,
          lineWidthPlus: 0
        }
      }
    }
  },
  exporting: {
    sourceWidth: 1200,
    scale: 1,
    chartOptions: {
      chart: {
        height: this && this.chartHeight
      },
      title: false
    }
  }
};

const SpectraViewer = ({ series }) => {
  const options = { ...SPECTRUM_CHARTOPTS, series };
  const chart = useRef();

  if (chart && chart.current) {
    window.chart = chart.current.chart;
  }

  console.log("HighchartsReact Render")
  return <HighchartsReact highcharts={Highcharts} options={options} ref={chart}/>;
};

export { SpectraViewer, SPECTRUM_CHARTOPTS };

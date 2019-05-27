import React, { useState, useEffect, useContext } from "react";
import Highcharts from "highcharts";
import {
  withHighcharts,
  HighchartsChart,
  Chart,
  Legend,
  XAxis,
  YAxis,
  AreaSplineSeries,
  Tooltip,
  Credits
} from "react-jsx-highcharts";
import applyExporting from "highcharts/modules/exporting";
import applyExportingData from "highcharts/modules/export-data";
import Switch from "@material-ui/core/Switch";
import FormGroup from "@material-ui/core/FormGroup";
import FormControlLabel from "@material-ui/core/FormControlLabel";
import { AppContext, initialize } from "./Store";
import { fetchSpectraList } from "./util";
import fixLogScale from "./fixLogScale";

applyExporting(Highcharts);
applyExportingData(Highcharts);
fixLogScale(Highcharts);

const FONTS =
  'Roboto, -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";';

const SPECTRUM_CHARTOPTS = {
  chart: {
    height: "350px",
    zoomType: "x",
    type: "areaspline",
    animation: { duration: 150 },
    panKey: "meta",
    panning: true
  },
  plotOptions: {
    areaspline: {},
    series: {
      fillOpacity: 0.45,
      animation: false,
      lineWidth: 0.4,
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
  credits: false,
  legend: {
    verticalAlign: "top",
    align: "right",
    x: -45,
    y: -1,
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
      fontFamily: FONTS,
      fontSize: "0.8rem"
    }
  },
  series: [],
  title: false,
  xAxis: {
    tickLength: 0,
    labels: {
      y: 15
    },
    min: null,
    max: null,
    crosshair: true
  },
  yAxis: {
    gridLineWidth: 0
    // title: false
    // reversed: true,
  },
  exporting: {
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
        menuItems: [
          "downloadPNG",
          "downloadPDF",
          "downloadSVG",
          "separator",
          "downloadCSV",
          "printChart",
          "separator",
          "reset"
          // "openInCloud"
          // "viewData"
        ]
      }
    }
  }
};

const SpectraViewer = () => {
  const [series, setSeries] = useState([]);
  const { state, dispatch } = useContext(AppContext);

  useEffect(() => {
    window.STATE = state;
    const updateSeriesData = async () => {
      const { currentSpectra } = state;
      let newSeries = series.filter(x => currentSpectra.includes(x.id));
      const seriesIDs = newSeries.map(x => x.id);
      const added = currentSpectra.filter(x => !seriesIDs.includes(x));

      if (added.length) {
        const newData = await fetchSpectraList(added);
        newSeries = [...newSeries, ...newData];
      }
      setSeries(newSeries);
    };
    updateSeriesData();
  }, [state.currentSpectra]); // eslint-disable-line react-hooks/exhaustive-deps

  const exportOptions = { ...SPECTRUM_CHARTOPTS.exporting };
  exportOptions.menuItemDefinitions = {
    reset: {
      onclick: () => {
        dispatch({ type: "RESET" });
        initialize(dispatch, false);
      },
      text: "Remove all spectra"
    }
  };

  const OD = num => (num <= 0 ? 10 : -Math.log10(num));
  const [logScale, setLogScale] = useState(false);
  const showOD = logScale || state.formState.F.filter(i => i.value).length > 0;

  return (
    <div>
      <HighchartsChart
        plotOptions={SPECTRUM_CHARTOPTS.plotOptions}
        navigation={SPECTRUM_CHARTOPTS.navigation}
        exporting={exportOptions}
      >
        <Chart {...SPECTRUM_CHARTOPTS.chart} />
        <Credits position={{ y: -45 }}>fpbase.org</Credits>
        <Legend {...SPECTRUM_CHARTOPTS.legend} />

        <XAxis {...SPECTRUM_CHARTOPTS.xAxis}>
          <XAxis.Title style={{ display: "none" }}>Wavelength</XAxis.Title>
        </XAxis>

        <YAxis
          {...SPECTRUM_CHARTOPTS.yAxis}
          reversed={logScale}
          min={0}
          max={logScale ? 6 : 1}
          labels={{ enabled: logScale }}
        >
          {series.map(serie => {
            let data = [...serie.data];
            if (logScale) {
              data = [...serie.data].map(([a, b]) => [a, OD(b)]);
            }
            return (
              <AreaSplineSeries
                key={serie.id}
                {...serie}
                data={data}
                threshold={logScale ? 10 : 0}
              />
            );
          })}
        </YAxis>

        <Tooltip {...SPECTRUM_CHARTOPTS.tooltip} />
      </HighchartsChart>

      {showOD && (
        <FormGroup
          style={{
            position: "absolute",
            top: 30,
            right: 0,
            opacity: logScale ? 0.8 : 0.5
          }}
        >
          <FormControlLabel
            // prettier-ignore
            control={(
              <Switch
                color="default"
                checked={logScale}
                onChange={(e,v)=> setLogScale(v)}
              />
)}
            label="OD"
          />
        </FormGroup>
      )}
    </div>
  );
};

export default withHighcharts(SpectraViewer, Highcharts);

import React, { useEffect, memo } from "react"
import { useQuery } from "@apollo/react-hooks"
import Highcharts from "highcharts"
import {
  withHighcharts,
  HighchartsChart,
  YAxis,
  Credits,
  Legend,
  Tooltip,
  provideAxis,
  Chart /* etc... */,
  XAxis,
} from "react-jsx-highcharts"
import applyExporting from "highcharts/modules/exporting"
import applyPatterns from "highcharts/modules/pattern-fill"
import applyExportingData from "highcharts/modules/export-data"
import addBoostModule from "highcharts/modules/boost"
import update from "immutability-helper"
import { css } from "@emotion/core"
import { BarLoader } from "react-spinners"
import gql from "graphql-tag"
import XRangePickers from "./XRangePickers"
import SpectrumSeries from "./SpectrumSeries"
import fixLogScale from "./fixLogScale"
import DEFAULT_OPTIONS from "./ChartOptions"
import NoData from "./NoData"
import useWindowWidth from "../useWindowWidth"
import useSpectraData from "../useSpectraData"

applyExporting(Highcharts)
applyExportingData(Highcharts)
applyPatterns(Highcharts)
addBoostModule(Highcharts)
fixLogScale(Highcharts)

const override = css`
  display: block;
  position: absolute;
  left: 40%;
  top: 50%;
  width: 20%;
  z-index: 10;
`

const calcHeight = width => {
  if (width < 600) return 275
  if (width < 960) return 325
  if (width < 1280) return 370
  if (width < 1920) return 400
  return 420
}

// export function spectraSorter(a, b) {
//   const SPECTRA_ORDER = ["P", "D", "F", "L", "C"]
//   if (SPECTRA_ORDER.indexOf(a.category) > SPECTRA_ORDER.indexOf(b.category))
//     return 1
//   return -1
// }

const {
  plotOptions: _plotOptions,
  xAxis: _xAxis,
  yAxis: _yAxis,
  chart: _chart,
  navigation: _navigation,
  exporting: _exporting,
  legend: _legend,
  tooltip: _tooltip,
  boost: _boost,
} = DEFAULT_OPTIONS

const BaseSpectraViewerContainer = React.memo(
  function BaseSpectraViewerContainer({
    ownerInfo,
    provideOptions,
    provideSpectra,
    provideOverlaps,
    provideHidden = [],
  }) {
    const { data } = useQuery(
      gql`
        query ChartOptions {
          chartOptions @client {
            showY
            showX
            showGrid
            areaFill
            logScale
            scaleEC
            scaleQY
            extremes
            shareTooltip
            palette
          }
          exNorm @client
        }
      `,
      { skip: provideOptions }
    )
    let chartOptions
    let normWave
    if (provideOptions) {
      chartOptions = provideOptions
      normWave = undefined
    } else {
      chartOptions = data.chartOptions
        ;[normWave] = data.exNorm
    }

    const yAxis = update(_yAxis, {
      labels: {
        enabled: { $set: chartOptions.showY || chartOptions.logScale },
      },
      gridLineWidth: { $set: chartOptions.showGrid ? 1 : 0 },
    })

    const xAxis = update(_xAxis, {
      labels: { enabled: { $set: chartOptions.showX } },
      gridLineWidth: { $set: chartOptions.showGrid ? 1 : 0 },
    })

    const tooltip = update(_tooltip, {
      shared: { $set: chartOptions.shareTooltip },
    })

    const spectraldata = useSpectraData(provideSpectra, provideOverlaps)
    return (
      <BaseSpectraViewer
        data={spectraldata}
        tooltip={tooltip}
        yAxis={yAxis}
        xAxis={xAxis}
        chartOptions={chartOptions}
        exNorm={+normWave}
        ownerInfo={ownerInfo}
        hidden={provideHidden}
      />
    )
  }
)

export const BaseSpectraViewer = memo(function BaseSpectraViewer({
  data,
  tooltip,
  yAxis,
  xAxis,
  exNorm,
  chartOptions,
  ownerInfo,
  hidden,
}) {
  const windowWidth = useWindowWidth()
  const numSpectra = data.length
  const owners = [...new Set(data.map(item => item.owner.slug))]
  const exData = data.filter(i => i.subtype === "EX" || i.subtype === "AB")
  const nonExData = data.filter(i => i.subtype !== "EX" && i.subtype !== "AB")

  let height = calcHeight(windowWidth) * (chartOptions.height || 1)
  let showPickers = numSpectra > 0 && !chartOptions.simpleMode
  if (chartOptions.zoomType !== undefined) {
    _chart.zoomType = chartOptions.zoomType
    showPickers = chartOptions.zoomType !== null
    // convert to no-op function
    xAxis.events.afterSetExtremes = () => { }
    xAxis.min = chartOptions.extremes[0]
    xAxis.max = chartOptions.extremes[1]
  }
  console.log(showPickers)
  const hChart = Highcharts.charts[0]
  let legendHeight
  if (hChart) {
    legendHeight = Highcharts.charts[0].legend.legendHeight || 0
    height += legendHeight
  }

  return (
    <div
      className="spectra-viewer"
      style={{ position: "relative", height: height }}
    >
      <span
        id="zoom-info"
        style={{
          display: "none",
          position: "absolute",
          fontWeight: 600,
          textAlign: "center",
          bottom: -1,
          width: "100%",
          zIndex: 10,
          fontSize: "0.7rem",
          color: "#bbb",
        }}
      />
      {numSpectra === 0 &&
        (chartOptions.simpleMode ? (
          <div className="sweet-loading">
            <BarLoader
              css={override}
              sizeUnit="px"
              size={100}
              color="#ccc"
              loading
            />
          </div>
        ) : (
          <NoData height={height} />
        ))}

      <ExNormNotice
        exNorm={exNorm}
        ownerInfo={ownerInfo}
        ecNorm={chartOptions.scaleEC}
        qyNorm={chartOptions.scaleQY}
      />
      <HighchartsChart
        plotOptions={_plotOptions}
        navigation={_navigation}
        exporting={_exporting}
        boost={_boost}
      >
        <Chart {..._chart} height={height} />
        <Legend {..._legend} />
        <Tooltip {...tooltip} />
        {/* the first Yaxis is for everything besides excitation data */}
        <YAxis
          id="yAx1"
          {...yAxis}
          reversed={chartOptions.logScale}
          max={chartOptions.logScale ? 6 : 1}
          min={0}
          gridLineWidth={chartOptions.showGrid && numSpectra > 0 ? 1 : 0}
          endOnTick={chartOptions.scaleEC}
          labels={{
            ...yAxis.labels,
            enabled: yAxis.labels.enabled && numSpectra > 0,
          }}
        >
          {nonExData.map(spectrum => (
            <SpectrumSeries
              exNorm={exNorm}
              spectrum={spectrum}
              key={spectrum.id}
              visible={!hidden.includes(spectrum.id)}
              ownerInfo={ownerInfo}
              ownerIndex={owners.indexOf(spectrum.owner.slug)}
              {...chartOptions}
            />
          ))}
        </YAxis>
        {/* a second axis for ex data, which may need to be scaled by EC */}
        <YAxis
          id="yAx2"
          {...yAxis}
          labels={{
            ...yAxis.labels,
            enabled: chartOptions.scaleEC,
            style: { fontWeight: 600, fontSize: "0.65rem" },
          }}
          opposite
          gridLineWidth={chartOptions.scaleEC && chartOptions.showGrid}
          maxPadding={0.0}
          reversed={chartOptions.logScale}
          max={chartOptions.scaleEC ? null : chartOptions.logScale ? 6 : 1}
          min={0}
          endOnTick={chartOptions.scaleEC}
        >
          {exData.length > 0 && chartOptions.scaleEC && (
            <YAxis.Title style={{ fontSize: "0.65rem" }}>
              Extinction Coefficient
            </YAxis.Title>
          )}
          {exData.map(spectrum => (
            <SpectrumSeries
              spectrum={spectrum}
              key={spectrum.id}
              ownerIndex={owners.indexOf(spectrum.owner.slug)}
              {...chartOptions}
            />
          ))}
        </YAxis>

        <XAxisWithRange
          options={xAxis}
          showPickers={showPickers}
        />
        <MyCredits
          axisId="yAx2"
          hide={numSpectra < 1 || chartOptions.simpleMode}
        />
      </HighchartsChart>
    </div>
  )
})

const MyCredits = provideAxis(function MyCredits({
  getAxis,
  getHighcharts,
  hide,
}) {
  useEffect(() => {
    const axis = getAxis()
    function shiftCredits() {
      const yShift = axis.object.chart.get("xAxis").axisTitleMargin
      axis.object.chart.credits.update({
        position: { y: -25 - yShift, x: -25 - axis.object.axisTitleMargin },
      })
    }
    const hChart = getHighcharts()
    hChart.addEvent(axis.object.chart, "redraw", shiftCredits)
    shiftCredits()
  }, []) // eslint-disable-line

  return (
    <Credits
      position={{ y: -45 }}
      href="https://www.fpbase.org/spectra"
      style={{ display: hide ? "none" : "block" }}
    >
      fpbase.org
    </Credits>
  )
})

export const XAxisWithRange = memo(function XAxisWithRange({
  options,
  showPickers,
}) {
  return (
    <>
      <XAxis {...options} lineWidth={showPickers ? 1 : 0} id="xAxis">
        <XAxis.Title style={{ display: "none" }}>Wavelength</XAxis.Title>
      </XAxis>
      {showPickers && (
        <XRangePickers
          axisId="xAxis"
          visible={showPickers && options.labels.enabled}
        />
      )}
    </>
  )
})

const ExNormNotice = memo(function ExNormNotice({
  exNorm,
  ecNorm,
  qyNorm,
  ownerInfo,
}) {
  const exNormed = ecNorm && Object.keys(ownerInfo).length > 0
  const emNormed = (exNorm || qyNorm) && Object.keys(ownerInfo).length > 0
  return (
    <div
      style={{
        position: "relative",
        top: -11,
        left: 20,
        zIndex: 1000,
        color: "rgba(200,0,0,0.45)",
        fontWeight: 600,
        fontSize: "0.82rem",
        height: 0,
      }}
    >
      {exNormed ? `EX NORMED TO EXT COEFF ${emNormed ? " & " : ""}` : ""}
      {emNormed && "EM NORMED TO "}
      {exNorm ? `${exNorm} EX${qyNorm ? " & " : ""}` : ""}
      {qyNorm && "QY"}
    </div>
  )
})

export const SpectraViewerContainer = withHighcharts(
  BaseSpectraViewerContainer,
  Highcharts
)
export const SpectraViewer = withHighcharts(BaseSpectraViewer, Highcharts)

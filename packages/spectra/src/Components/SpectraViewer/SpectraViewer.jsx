import Highcharts from "highcharts"
import React, { memo, useEffect } from "react"
import {
  Chart /* etc... */,
  Credits,
  HighchartsChart,
  HighchartsProvider,
  Legend,
  Tooltip,
  useAxis,
  useHighcharts,
  XAxis,
  YAxis,
} from "react-jsx-highcharts"
import "highcharts/modules/exporting"
import "highcharts/modules/offline-exporting"
import "highcharts/modules/pattern-fill"
import "highcharts/modules/export-data"
import "highcharts/modules/boost"
import LinearProgress from "@mui/material/LinearProgress"
import { defaultChartOptions } from "../../defaults"
import useSpectraData from "../../hooks/useSpectraData"
import { useSpectraStore } from "../../store/spectraStore"
import useWindowWidth from "../useWindowWidth"
import DEFAULT_OPTIONS from "./ChartOptions"
import fixLogScale from "./fixLogScale"
import NoData from "./NoData"
import SpectrumSeries from "./SpectrumSeries"
import { XAxisRangeInputs } from "./XAxisRangeInputs"

fixLogScale(Highcharts)

const calcHeight = (width) => {
  if (width < 600) return 275
  if (width < 960) return 325
  if (width < 1280) return 370
  if (width < 1920) return 400
  return 420
}

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

const BaseSpectraViewerContainer = React.memo(function BaseSpectraViewerContainer({
  ownerInfo,
  provideState,
}) {
  // Get state from Zustand store or use provided state
  const storeActiveSpectra = useSpectraStore((state) => state.activeSpectra)
  const storeActiveOverlaps = useSpectraStore((state) => state.activeOverlaps)
  const storeHiddenSpectra = useSpectraStore((state) => state.hiddenSpectra)
  const storeChartOptions = useSpectraStore((state) => state.chartOptions)
  const storeExNorm = useSpectraStore((state) => state.exNorm)

  // Use provided state or fall back to store
  const activeSpectra = provideState?.activeSpectra ?? storeActiveSpectra
  const activeOverlaps = provideState?.activeOverlaps ?? storeActiveOverlaps
  const hiddenSpectra = provideState?.hiddenSpectra ?? storeHiddenSpectra

  // Merge provided chartOptions with defaults to ensure all required fields are present
  const chartOptions = provideState?.chartOptions
    ? { ...defaultChartOptions, ...provideState.chartOptions }
    : storeChartOptions

  // Always call useSpectraData hook before any returns (Rules of Hooks)
  const spectraldata = useSpectraData(activeSpectra, activeOverlaps)

  // Safely extract normWave from exNorm array, handling non-array values
  // exNorm should be [normWave, normID] but may be corrupted/malformed
  let normWave
  if (provideState?.chartOptions) {
    // If state is provided, don't use exNorm from store
    normWave = undefined
  } else if (Array.isArray(storeExNorm)) {
    ;[normWave] = storeExNorm
  } else {
    normWave = null
  }

  const yAxis = {
    ..._yAxis,
    labels: {
      ..._yAxis.labels,
      enabled: chartOptions.showY || chartOptions.logScale,
    },
    gridLineWidth: chartOptions.showGrid ? 1 : 0,
  }

  const updateChartOptions = useSpectraStore((state) => state.updateChartOptions)

  const xAxis = {
    ..._xAxis,
    labels: {
      ..._xAxis.labels,
      enabled: chartOptions.showX && activeSpectra.length > 0,
    },
    gridLineWidth: chartOptions.showGrid ? 1 : 0,
    events: {
      ..._xAxis.events,
      afterSetExtremes: (event) => {
        // Call original handler for zoom-info display
        if (_xAxis.events?.afterSetExtremes) {
          _xAxis.events.afterSetExtremes(event)
        }

        const { min, max, userMin, userMax, dataMin, dataMax, trigger } = event

        // Handle reset case: both min and max are null OR extremes match full data range
        if (min === null && max === null) {
          if (chartOptions.extremes !== null) {
            updateChartOptions({ extremes: null })
          }
          return
        }

        // If extremes match the full data range (with tolerance for Highcharts padding), treat as "no zoom" (reset)
        // This handles the case where Highcharts reset button sets extremes to data range with padding
        const dataRange = dataMax - dataMin
        const tolerance = dataRange * 0.02 // 2% tolerance for padding
        const isFullDataRange =
          Math.abs(min - dataMin) <= tolerance && Math.abs(max - dataMax) <= tolerance
        if (isFullDataRange) {
          if (chartOptions.extremes !== null) {
            updateChartOptions({ extremes: null })
          }
          return
        }

        // Only save extremes if this was a user-initiated zoom
        // Skip auto-fitting and programmatic updates (our own XAxisRangeInputs component)
        // event.trigger can be: 'zoom', 'navigator', 'rangeSelectorButton', 'rangeSelectorInput', undefined
        const isUserZoom = trigger === "zoom"
        if (!isUserZoom) {
          return
        }

        // Zoom case - use userMin/userMax (actual zoom values set by user)
        const newMin = userMin ?? min
        const newMax = userMax ?? max

        // Only update if extremes actually changed (prevents infinite loops)
        const currentMin = chartOptions.extremes?.[0]
        const currentMax = chartOptions.extremes?.[1]

        if (currentMin !== newMin || currentMax !== newMax) {
          updateChartOptions({
            extremes: [newMin ?? null, newMax ?? null],
          })
        }
      },
    },
  }

  const tooltip = {
    ..._tooltip,
    shared: chartOptions.shareTooltip,
  }

  return (
    <BaseSpectraViewer
      data={spectraldata}
      tooltip={tooltip}
      yAxis={yAxis}
      xAxis={xAxis}
      chartOptions={chartOptions}
      exNorm={+normWave}
      ownerInfo={ownerInfo}
      hidden={hiddenSpectra}
    />
  )
})

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
  const owners = [
    ...new Set(data.map((item) => item.owner?.slug).filter((slug) => slug !== undefined)),
  ]
  const exData = data.filter((i) => i.subtype === "EX" || i.subtype === "AB")
  const nonExData = data.filter((i) => i.subtype !== "EX" && i.subtype !== "AB")

  const height = calcHeight(windowWidth) * (chartOptions.height || 1)
  if (chartOptions.zoomType !== undefined) {
    _chart.zoomType = chartOptions.zoomType
    // convert to no-op function
    xAxis.events.afterSetExtremes = () => {}
    // Handle new extremes format: [number | null, number | null] | null
    if (chartOptions.extremes) {
      xAxis.min = chartOptions.extremes[0] ?? undefined
      xAxis.max = chartOptions.extremes[1] ?? undefined
    }
  }
  // Note: legendHeight is already accounted for by Highcharts internally
  // Adding it here causes reflow issues when toggling series visibility

  return (
    <div
      id="spectra-viewer-container"
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
            <LinearProgress
              sx={{
                position: "absolute",
                left: "40%",
                top: "50%",
                width: "20%",
                height: 4,
                zIndex: 10,
              }}
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
      <HighchartsProvider Highcharts={Highcharts}>
        <HighchartsChart
          plotOptions={_plotOptions}
          navigation={_navigation}
          exporting={_exporting}
          lang={{ noData: "" }}
          accessibility={{ enabled: false }}
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
            {nonExData.map((spectrum) => (
              <SpectrumSeries
                exNorm={exNorm}
                spectrum={spectrum}
                key={spectrum.id}
                visible={!hidden.includes(spectrum.id)}
                ownerInfo={ownerInfo}
                ownerIndex={spectrum.owner?.slug ? owners.indexOf(spectrum.owner.slug) : -1}
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
              <YAxis.Title style={{ fontSize: "0.65rem" }}>Extinction Coefficient</YAxis.Title>
            )}
            {exData.map((spectrum) => (
              <SpectrumSeries
                spectrum={spectrum}
                key={spectrum.id}
                visible={!hidden.includes(spectrum.id)}
                ownerIndex={spectrum.owner?.slug ? owners.indexOf(spectrum.owner.slug) : -1}
                {...chartOptions}
              />
            ))}
            <MyCredits hide={numSpectra < 1 || chartOptions.simpleMode} />
          </YAxis>

          <XAxis {...xAxis} lineWidth={numSpectra > 0 ? 1 : 0} id="xAxis">
            <XAxis.Title style={{ display: "none" }}>Wavelength</XAxis.Title>
            <XAxisRangeInputs enabled={chartOptions.showX && numSpectra > 0} />
          </XAxis>
        </HighchartsChart>
      </HighchartsProvider>
    </div>
  )
})

const MyCredits = function MyCredits({ hide }) {
  const axis = useAxis()
  const Highcharts = useHighcharts()

  useEffect(() => {
    if (!axis || !axis.object || !Highcharts) return

    function shiftCredits() {
      const yShift = axis.object.chart.get("xAxis").axisTitleMargin
      axis.object.chart.credits.update({
        position: { y: -25 - yShift, x: -25 - axis.object.axisTitleMargin },
      })
    }
    Highcharts.addEvent(axis.object.chart, "redraw", shiftCredits)
    shiftCredits()
  }, [axis, Highcharts])

  return (
    <Credits
      position={{ y: -45 }}
      href="https://www.fpbase.org/spectra"
      style={{ display: hide ? "none" : "block" }}
    >
      fpbase.org
    </Credits>
  )
}

export const XAxisWithRange = memo(function XAxisWithRange({ options }) {
  return (
    <XAxis {...options} id="xAxis">
      <XAxis.Title style={{ display: "none" }}>Wavelength</XAxis.Title>
    </XAxis>
  )
})

const ExNormNotice = memo(function ExNormNotice({ exNorm, ecNorm, qyNorm, ownerInfo }) {
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

export const SpectraViewerContainer = BaseSpectraViewerContainer
export const SpectraViewer = BaseSpectraViewer

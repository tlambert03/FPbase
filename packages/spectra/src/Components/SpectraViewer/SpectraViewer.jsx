import Highcharts from "highcharts"
import { HighchartsReact } from "highcharts-react-official"
import React, { memo, useEffect, useMemo, useRef, useState } from "react"
import "highcharts/modules/exporting"
import "highcharts/modules/offline-exporting"
import "highcharts/modules/pattern-fill"
import "highcharts/modules/export-data"
import LinearProgress from "@mui/material/LinearProgress"
import { defaultChartOptions } from "../../defaults"
import useSpectraData from "../../hooks/useSpectraData"
import { useSpectraBatch } from "../../hooks/useSpectraQueries"
import { useSpectraStore } from "../../store/spectraStore"
import useWindowWidth from "../useWindowWidth"
import { AxisContext } from "./axisContext"
import DEFAULT_OPTIONS from "./ChartOptions"
import fixLogScale from "./fixLogScale"
import NoData from "./NoData"
import { exNormSpectrumId, spectrumSeriesOptions } from "./SpectrumSeries"
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

  const updateChartOptions = useSpectraStore((state) => state.updateChartOptions)

  // With provided state (protein pages, embeds) the viewer is isolated from the persisted
  // store: zoom extremes live in local state instead of being written to the store
  const isolated = provideState != null
  const [localExtremes, setLocalExtremes] = useState(provideState?.chartOptions?.extremes ?? null)
  const setExtremes = isolated ? setLocalExtremes : (extremes) => updateChartOptions({ extremes })

  // Merge provided chartOptions with defaults to ensure all required fields are present
  const chartOptions = isolated
    ? { ...defaultChartOptions, ...provideState.chartOptions, extremes: localExtremes }
    : storeChartOptions

  // Always call useSpectraData hook before any returns (Rules of Hooks)
  const spectraldata = useSpectraData(activeSpectra, activeOverlaps)

  // Safely extract normWave from exNorm array, handling non-array values
  // exNorm should be [normWave, normID] but may be corrupted/malformed
  let normWave
  if (isolated) {
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

        // Only process user-initiated zoom events (not programmatic updates)
        // event.trigger can be: 'zoom', 'navigator', 'rangeSelectorButton', 'rangeSelectorInput', undefined
        // This check MUST come first to prevent programmatic updates (like restoring from sessionStorage)
        // from incorrectly clearing extremes
        const isUserZoom = trigger === "zoom"
        if (!isUserZoom) {
          return
        }

        // Handle reset case: both min and max are null
        if (min === null && max === null) {
          if (chartOptions.extremes !== null) {
            setExtremes(null)
          }
          return
        }

        // If extremes match the full data range (with tolerance for Highcharts padding), treat as "no zoom" (reset)
        // This handles the case where user zooms out to full range or clicks reset zoom button
        const dataRange = dataMax - dataMin
        const tolerance = dataRange * 0.02 // 2% tolerance for padding
        const isFullDataRange =
          Math.abs(min - dataMin) <= tolerance && Math.abs(max - dataMax) <= tolerance
        if (isFullDataRange) {
          if (chartOptions.extremes !== null) {
            setExtremes(null)
          }
          return
        }

        // Zoom case - use userMin/userMax (actual zoom values set by user)
        const newMin = userMin ?? min
        const newMax = userMax ?? max

        // Only update if extremes actually changed (prevents infinite loops)
        const currentMin = chartOptions.extremes?.[0]
        const currentMax = chartOptions.extremes?.[1]

        if (currentMin !== newMin || currentMax !== newMax) {
          setExtremes([newMin ?? null, newMax ?? null])
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
      onExtremesChange={isolated ? setLocalExtremes : undefined}
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
  onExtremesChange,
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

  const exNormIds = useMemo(
    () => [...new Set(data.map((s) => exNormSpectrumId(s, ownerInfo, exNorm)).filter(Boolean))],
    [data, ownerInfo, exNorm]
  )
  const { data: exNormSpectra } = useSpectraBatch(exNormIds)

  const series = useMemo(() => {
    const exById = Object.fromEntries(exNormSpectra.map((s) => [String(s.id), s]))
    const toSeries = (spectrum, yAxis, withExNorm) =>
      spectrumSeriesOptions(
        spectrum,
        { ...chartOptions, exNorm: withExNorm ? exNorm : undefined },
        {
          exSpectrum: withExNorm ? exById[exNormSpectrumId(spectrum, ownerInfo, exNorm)] : null,
          ownerIndex: spectrum.owner?.slug ? owners.indexOf(spectrum.owner.slug) : -1,
          visible: !hidden.includes(spectrum.id),
          yAxis,
        }
      )
    const valid = (s) => s?.id && s.data
    return [
      ...nonExData.filter(valid).map((s) => toSeries(s, "yAx1", true)),
      ...exData.filter(valid).map((s) => toSeries(s, "yAx2", false)),
    ]
  }, [data, chartOptions, exNorm, ownerInfo, hidden, exNormSpectra])

  const hideCredits = numSpectra < 1 || chartOptions.simpleMode
  const options = {
    chart: { ..._chart, height },
    title: { text: null },
    subtitle: { text: null },
    plotOptions: _plotOptions,
    navigation: _navigation,
    exporting: _exporting,
    lang: { noData: "" },
    accessibility: { enabled: false },
    legend: { ..._legend, enabled: _legend.enabled ?? true },
    tooltip: { ...tooltip, enabled: tooltip.enabled ?? true },
    credits: {
      enabled: true,
      text: "fpbase.org",
      href: "https://www.fpbase.org/spectra",
      position: { y: -45 },
      style: { display: hideCredits ? "none" : "block" },
    },
    yAxis: [
      // the first Yaxis is for everything besides excitation data
      {
        id: "yAx1",
        title: { text: null },
        ...yAxis,
        reversed: chartOptions.logScale,
        max: chartOptions.logScale ? 6 : 1,
        min: 0,
        gridLineWidth: chartOptions.showGrid && numSpectra > 0 ? 1 : 0,
        endOnTick: chartOptions.scaleEC,
        labels: { ...yAxis.labels, enabled: yAxis.labels.enabled && numSpectra > 0 },
      },
      // a second axis for ex data, which may need to be scaled by EC
      {
        id: "yAx2",
        ...yAxis,
        title: {
          text: exData.length > 0 && chartOptions.scaleEC ? "Extinction Coefficient" : null,
          style: { fontSize: "0.65rem" },
        },
        labels: {
          ...yAxis.labels,
          enabled: chartOptions.scaleEC,
          style: { fontWeight: 600, fontSize: "0.65rem" },
        },
        opposite: true,
        gridLineWidth: chartOptions.scaleEC && chartOptions.showGrid,
        maxPadding: 0.0,
        reversed: chartOptions.logScale,
        max: chartOptions.scaleEC ? null : chartOptions.logScale ? 6 : 1,
        min: 0,
        endOnTick: chartOptions.scaleEC,
      },
    ],
    xAxis: [
      {
        id: "xAxis",
        ...xAxis,
        title: { text: "Wavelength", style: { display: "none" } },
        lineWidth: numSpectra > 0 ? 1 : 0,
      },
    ],
    series,
  }

  // read the chart from the ref, not `callback`: exporting builds temporary chart copies
  // and runs the original chart's callback for them too
  const chartRef = useRef(null)
  const [chart, setChart] = useState(null)
  useEffect(() => {
    setChart(chartRef.current?.chart ?? null)
  }, [])
  const xAxisContext = useMemo(
    () => (chart ? { object: chart.get("xAxis"), id: "xAxis" } : null),
    [chart]
  )

  // keep the credits clear of the axis titles
  useEffect(() => {
    if (!chart) return
    const shiftCredits = () => {
      const yShift = chart.get("xAxis").axisTitleMargin
      chart.credits?.update({
        position: { y: -25 - yShift, x: -25 - chart.get("yAx2").axisTitleMargin },
      })
    }
    shiftCredits()
    return Highcharts.addEvent(chart, "redraw", shiftCredits)
  }, [chart])

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
      <HighchartsReact
        highcharts={Highcharts}
        options={options}
        ref={chartRef}
        containerProps={{ className: "chart" }}
      />
      {xAxisContext && (
        <AxisContext.Provider value={xAxisContext}>
          <XAxisRangeInputs
            enabled={chartOptions.showX && numSpectra > 0}
            extremes={chartOptions.extremes}
            onExtremesChange={onExtremesChange}
          />
        </AxisContext.Provider>
      )}
    </div>
  )
})

const ExNormNotice = memo(function ExNormNotice({ exNorm, ecNorm, qyNorm, ownerInfo = {} }) {
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

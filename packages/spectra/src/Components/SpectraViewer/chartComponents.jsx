/**
 * Minimal React bindings for a single Highcharts chart.
 *
 * Mirrors the subset of react-jsx-highcharts that the spectra viewer used: each
 * component adds/updates/removes its own piece of one chart instance, and redraws
 * are batched to one per animation frame.
 */
import Highcharts from "highcharts"
import {
  createContext,
  memo,
  useContext,
  useEffect,
  useLayoutEffect,
  useRef,
  useState,
} from "react"

const ChartContext = createContext(null)
/** @type {import("react").Context<{ object: import("highcharts").Axis, id: string } | null>} */
const AxisContext = createContext(null)

let nextId = 0
const uniqueId = (prefix) => `${prefix}-${++nextId}`

// props whose values changed since the last render (Object.is), or false if none
function useModifiedProps(props) {
  const prev = useRef(null)
  useEffect(() => {
    prev.current = props
  })
  const modified = {}
  for (const key of Object.keys(props)) {
    if (!prev.current || !Object.is(props[key], prev.current[key])) modified[key] = props[key]
  }
  return Object.keys(modified).length > 0 ? modified : false
}

function rafDebounce(fn) {
  let queued
  return () => {
    if (queued) cancelAnimationFrame(queued)
    queued = requestAnimationFrame(fn)
  }
}

export function HighchartsChart({ children = null, plotOptions, ...options }) {
  const domRef = useRef(null)
  const [provided, setProvided] = useState(null)

  useLayoutEffect(() => {
    const chart = Highcharts.chart(domRef.current, {
      title: { text: null },
      subtitle: { text: null },
      legend: { enabled: false },
      tooltip: { enabled: false },
      credits: { enabled: false },
      series: [],
      xAxis: [],
      yAxis: [],
      plotOptions,
      ...options,
    })
    const needsRedraw = rafDebounce(() => {
      if (!chart.__destroyed) chart.redraw()
    })
    setProvided({ object: chart, needsRedraw })
  }, [])

  useEffect(
    () => () => {
      if (provided) {
        provided.object.__destroyed = true
        provided.object.destroy()
      }
    },
    [provided]
  )

  const prevPlotOptions = useRef(plotOptions)
  useEffect(() => {
    if (provided && !Object.is(prevPlotOptions.current, plotOptions)) {
      provided.object.update({ plotOptions }, false)
      provided.needsRedraw()
    }
    prevPlotOptions.current = plotOptions
  })

  return (
    <div className="chart" ref={domRef}>
      {provided && <ChartContext.Provider value={provided}>{children}</ChartContext.Provider>}
    </div>
  )
}

const useChart = () => useContext(ChartContext)

export const Chart = memo(function Chart({ width, height, ...options }) {
  const chart = useChart()
  const mounted = useRef(false)
  const modified = useModifiedProps(options)

  useEffect(() => {
    if (!(width === undefined && height === undefined)) chart.object.setSize(width, height)
  }, [width, height])

  useEffect(() => {
    if (!mounted.current) {
      chart.object.update({ chart: options }, false)
      chart.needsRedraw()
      mounted.current = true
    } else if (modified !== false) {
      chart.object.update({ chart: modified }, false)
      chart.needsRedraw()
    }
  })
  return null
})

// Legend / Credits: push modified options; disable on unmount
function useChartSection(props, update) {
  const chart = useChart()
  const modified = useModifiedProps(props)
  useEffect(() => {
    if (modified !== false) {
      update(chart.object, modified)
      chart.needsRedraw()
    }
  })
  useEffect(
    () => () => {
      if (!chart.object.__destroyed) update(chart.object, { enabled: false })
      chart.needsRedraw()
    },
    []
  )
}

export const Legend = memo(function Legend({ enabled = true, ...options }) {
  useChartSection({ enabled, ...options }, (chart, config) =>
    chart.update({ legend: config }, false)
  )
  return null
})

export function Credits({ enabled = true, children, ...options }) {
  useChartSection({ enabled, text: children, ...options }, (chart, config) => {
    const { text, ...rest } = config
    chart.addCredits(text ? config : rest, true)
  })
  return null
}

export const Tooltip = memo(function Tooltip({ enabled = true, ...options }) {
  const chart = useChart()
  const props = { enabled, ...options }
  const modified = useModifiedProps(props)
  const mounted = useRef(false)
  useEffect(() => {
    if (!mounted.current) {
      chart.object.update({ tooltip: { ...Highcharts.defaultOptions?.tooltip, ...props } })
      mounted.current = true
    } else if (modified !== false) {
      chart.object.update({ tooltip: modified })
    }
  })
  useEffect(
    () => () => {
      if (!chart.object.__destroyed) chart.object.update({ tooltip: { enabled: false } })
    },
    []
  )
  return null
})

function Axis({ isX, id, children = null, ...options }) {
  const chart = useChart()
  const [axis, setAxis] = useState(null)
  const modified = useModifiedProps(options)

  useEffect(() => {
    const created = chart.object.addAxis(
      { id: id ?? uniqueId("axis"), title: { text: null }, ...options },
      isX,
      false
    )
    setAxis({ object: created, id: created.options.id })
    chart.needsRedraw()
    return () => {
      try {
        created.remove(false)
      } catch {
        // already removed along with the chart
      }
      chart.needsRedraw()
    }
  }, [])

  useEffect(() => {
    if (axis && modified !== false) {
      axis.object.update(modified, false)
      chart.needsRedraw()
    }
  })

  if (!axis) return null
  return <AxisContext.Provider value={axis}>{children}</AxisContext.Provider>
}

export const useAxis = () => useContext(AxisContext)

const AxisTitle = memo(function AxisTitle({ children: text, ...options }) {
  const axis = useAxis()
  useEffect(() => {
    if (axis) axis.object.setTitle({ text, ...options }, true)
  })
  useEffect(
    () => () => {
      try {
        axis?.object.setTitle({ text: null }, true)
      } catch {
        // axis already removed
      }
    },
    [axis]
  )
  return null
})

export const XAxis = ({ type = "linear", ...props }) => <Axis type={type} {...props} isX />
XAxis.Title = AxisTitle
export const YAxis = ({ type = "linear", ...props }) => <Axis type={type} {...props} isX={false} />
YAxis.Title = AxisTitle

const EMPTY = []

export const Series = memo(function Series({ data = EMPTY, visible = true, ...options }) {
  const chart = useChart()
  const axis = useAxis()
  const seriesRef = useRef(null)
  const prev = useRef(null)

  useEffect(() => {
    if (!axis) return
    const series = chart.object.addSeries(
      { id: uniqueId("series"), data, visible, ...options, [axis.object.coll]: axis.id },
      false
    )
    seriesRef.current = series
    chart.needsRedraw()
    return () => {
      try {
        series.remove(false)
      } catch {
        // already removed along with its axis
      }
      seriesRef.current = null
      chart.needsRedraw()
    }
  }, [axis])

  useEffect(() => {
    const series = seriesRef.current
    const last = prev.current
    prev.current = { data, visible, options }
    if (!series || !last) return
    let redraw = false
    if (!Object.is(data, last.data)) {
      series.setData(data, false)
      redraw = true
    }
    if (visible !== last.visible) {
      series.setVisible(visible, false)
      redraw = true
    }
    const modified = {}
    for (const key of Object.keys(options)) {
      if (!Object.is(options[key], last.options[key])) modified[key] = options[key]
    }
    if (Object.keys(modified).length > 0) {
      series.update(modified, false)
      redraw = true
    }
    if (redraw) chart.needsRedraw()
  })
  return null
})

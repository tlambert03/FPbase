import { Tooltip } from "@mui/material"
import Input from "@mui/material/Input"
import { withStyles } from "@mui/styles"
import { useEffect, useRef, useState } from "react"
import { useAxis, useHighcharts } from "react-jsx-highcharts"
import { useSpectraStore } from "../../store/spectraStore"

const LightTooltip = withStyles((theme) => ({
  tooltip: {
    backgroundColor: theme.palette.common.white,
    color: "rgba(0, 0, 0, 0.87)",
    boxShadow: theme.shadows[1],
    fontSize: 12,
    margin: "0 13px 7px",
  },
}))(Tooltip)

const CLASSES = {
  minInput: {
    fontWeight: "bold",
    fontSize: "0.75rem",
    width: 30,
    position: "absolute",
  },
  maxInput: {
    position: "absolute",
    fontWeight: "bold",
    fontSize: "0.75rem",
    width: 30,
  },
}

const XRangePickers = ({ visible }) => {
  const axis = useAxis()
  const Highcharts = useHighcharts()
  const storeExtremes = useSpectraStore((state) => state.chartOptions.extremes)
  const updateChartOptions = useSpectraStore((state) => state.updateChartOptions)

  // Local state for input values (allows typing without immediate updates)
  const [minValue, setMinValue] = useState("")
  const [maxValue, setMaxValue] = useState("")

  // Refs for positioning
  const minNode = useRef()
  const maxNode = useRef()

  // Initialize from store extremes on mount or when store extremes change
  useEffect(() => {
    if (!storeExtremes) return
    const [min, max] = storeExtremes
    if (min || max) {
      setMinValue(min || "")
      setMaxValue(max || "")

      // Apply extremes to axis and show reset zoom button
      // This handles the case when extremes are restored from sessionStorage on page load
      if (axis?.object && min && max) {
        axis.setExtremes(min, max)
        axis.object.chart.showResetZoom()
      }
    }
  }, [storeExtremes, axis])

  // Sync local state with axis extremes when axis changes (from zooming, etc)
  useEffect(() => {
    if (!axis?.object) return

    const handleExtremesChange = () => {
      const extremes = axis.object.getExtremes()
      // Store user-set extremes (null if autoscaled)
      const userMin = extremes.userMin ? Math.round(extremes.min) : null
      const userMax = extremes.userMax ? Math.round(extremes.max) : null

      // Display current extremes (whether user-set or autoscaled)
      const displayMin = Math.round(extremes.min)
      const displayMax = Math.round(extremes.max)

      setMinValue(String(displayMin))
      setMaxValue(String(displayMax))

      // Update store for persistence
      // Normalize [null, null] to just null (no zoom = null, not [null, null])
      const newExtremes = userMin === null && userMax === null ? null : [userMin, userMax]
      updateChartOptions({ extremes: newExtremes })
    }

    const handleRedraw = () => {
      // Position inputs relative to y-axis margins
      if (!axis.object.labelGroup || !minNode.current) return

      let leftPad = -5
      const yAx1 = axis.object.chart.get("yAx1")
      if (yAx1) leftPad += yAx1.axisTitleMargin

      let rightPad = 0
      const yAx2 = axis.object.chart.get("yAx2")
      if (yAx2) rightPad += yAx2.axisTitleMargin

      minNode.current.parentElement.style.left = `${leftPad}px`
      maxNode.current.parentElement.style.right = `${rightPad}px`

      // Hide axis labels that are too close to our inputs
      const extremes = axis.object.getExtremes()
      const threshold = 0.43 * axis.object.tickInterval

      axis.object.labelGroup.element.childNodes.forEach((node) => {
        const value = Number.parseFloat(node.textContent)
        const minDist = Math.abs(value - extremes.min)
        const maxDist = Math.abs(value - extremes.max)
        node.style.display = Math.min(minDist, maxDist) < threshold ? "none" : "block"
      })
    }

    if (!Highcharts) return

    Highcharts.addEvent(axis.object, "afterSetExtremes", handleExtremesChange)
    Highcharts.addEvent(axis.object.chart, "redraw", handleRedraw)

    // Initialize values
    handleExtremesChange()
    handleRedraw()

    return () => {
      try {
        Highcharts.removeEvent(axis.object, "afterSetExtremes", handleExtremesChange)
        Highcharts.removeEvent(axis.object.chart, "redraw", handleRedraw)
      } catch {
        // Axis/chart may be destroyed
      }
    }
  }, [axis, Highcharts, updateChartOptions])

  // Apply extremes to the axis
  const applyExtremes = (newMin, newMax) => {
    if (!axis?.object) return

    const min = newMin ? Number(newMin) : null
    const max = newMax ? Number(newMax) : null

    axis.setExtremes(min, max)

    // Show reset zoom button if both are set
    if (min && max) {
      axis.object.chart.showResetZoom()
    }
  }

  const handleMinChange = (e) => {
    setMinValue(e.target.value)
  }

  const handleMaxChange = (e) => {
    setMaxValue(e.target.value)
  }

  const handleMinBlur = () => {
    applyExtremes(minValue, maxValue)
  }

  const handleMaxBlur = () => {
    applyExtremes(minValue, maxValue)
  }

  const handleKeyDown = (e) => {
    if (e.key === "Enter") {
      applyExtremes(minValue, maxValue)
      e.target.select()
    }
  }

  if (!axis?.object) return null

  const axisExtremes = axis.getExtremes()

  // Color logic: gray if autoscaled, red if cropping data, blue if zoomed but showing all data
  const minColor = !axisExtremes.userMin
    ? "#444"
    : axisExtremes.dataMin < axisExtremes.min
      ? "#B1191E"
      : "#5F67CE"

  const maxColor = !axisExtremes.userMax
    ? "#444"
    : axisExtremes.dataMax > axisExtremes.max
      ? "#B1191E"
      : "#5F67CE"

  return (
    <div
      className="x-range-pickers"
      style={{
        height: 0,
        position: "relative",
        bottom: 38,
        display: visible ? "block" : "none",
      }}
    >
      <LightTooltip
        title="Type to change min, clear to autoscale"
        placement="top-end"
        TransitionProps={{ timeout: { enter: 150, exit: 400 } }}
      >
        <Input
          name="min"
          type="text"
          placeholder={`${Math.round(axisExtremes.dataMin) || ""}`}
          value={minValue}
          inputRef={minNode}
          onChange={handleMinChange}
          onKeyDown={handleKeyDown}
          onBlur={handleMinBlur}
          style={{ ...CLASSES.minInput, color: minColor }}
          inputProps={{ style: { textAlign: "center" } }}
        />
      </LightTooltip>
      <LightTooltip
        title="Type to change max, clear to autoscale"
        placement="top-start"
        TransitionProps={{ timeout: { enter: 150, exit: 400 } }}
      >
        <Input
          name="max"
          type="text"
          placeholder={`${Math.round(axisExtremes.dataMax) || ""}`}
          value={maxValue}
          inputRef={maxNode}
          onChange={handleMaxChange}
          onKeyDown={handleKeyDown}
          onBlur={handleMaxBlur}
          style={{ ...CLASSES.maxInput, color: maxColor }}
          inputProps={{ style: { textAlign: "center" } }}
        />
      </LightTooltip>
    </div>
  )
}

export default XRangePickers

import { Tooltip } from "@mui/material"
import Input from "@mui/material/Input"
import { withStyles } from "@mui/styles"
import React, { useEffect, useRef } from "react"
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
    color: "#444",
    position: "absolute",
  },
  maxInput: {
    position: "absolute",
    fontWeight: "bold",
    fontSize: "0.75rem",
    width: 30,
    color: "#444",
  },
}

const XRangePickers = ({ visible }) => {
  const axis = useAxis()
  const Highcharts = useHighcharts()
  const extremes = useSpectraStore((state) => state.chartOptions.extremes)
  const updateChartOptions = useSpectraStore((state) => state.updateChartOptions)
  const [min, max] = extremes || [null, null]
  const minNode = useRef()
  const maxNode = useRef()

  useEffect(() => {
    if (!axis || !axis.object) return
    if (min || max) {
      // Only call setExtremes if the values are actually different from current axis extremes
      // This prevents infinite loop with handleAfterSetExtremes
      const current = axis.object.getExtremes()
      const targetMin = min && Math.round(min)
      const targetMax = max && Math.round(max)
      const currentMin = current.userMin && Math.round(current.min)
      const currentMax = current.userMax && Math.round(current.max)

      // Only update if the extremes have actually changed
      if (targetMin !== currentMin || targetMax !== currentMax) {
        axis.setExtremes(targetMin, targetMax)
        if (min && max) {
          axis.object.chart.showResetZoom()
        }
      }
    }
  }, [axis, min, max])

  useEffect(() => {
    if (!axis || !axis.object || !Highcharts) return

    function handleAfterSetExtremes() {
      const e = axis.object.getExtremes()
      if (e) {
        const newExtremes = [e.userMin && Math.round(e.min), e.userMax && Math.round(e.max)]
        // Only update store if extremes actually changed (prevents infinite loop)
        const currentExtremes = useSpectraStore.getState().chartOptions.extremes || [null, null]
        if (newExtremes[0] !== currentExtremes[0] || newExtremes[1] !== currentExtremes[1]) {
          updateChartOptions({ extremes: newExtremes })
        }
      }
    }

    function positionInputs() {
      if (axis.object.labelGroup && minNode.current) {
        let leftPad = -5
        if (axis.object.chart.get("yAx1")) {
          leftPad += +axis.object.chart.get("yAx1").axisTitleMargin
        }
        let rightPad = 0
        if (axis.object.chart.get("yAx2")) {
          rightPad += +axis.object.chart.get("yAx2").axisTitleMargin
        }
        minNode.current.parentElement.style.left = `${leftPad}px`
        maxNode.current.parentElement.style.right = `${rightPad}px`
        axis.object.labelGroup.element.childNodes.forEach((node) => {
          node.style.display = "block"
        })
        const { min: exMin, max: exMax } = axis.getExtremes()
        axis.object.labelGroup.element.childNodes.forEach((node) => {
          if (
            Math.min(Math.abs(node.textContent - exMin), Math.abs(node.textContent - exMax)) <
            0.43 * axis.object.tickInterval
          ) {
            node.style.display = "none"
          }
        })
      }
    }

    Highcharts.addEvent(axis.object.chart, "redraw", positionInputs)
    Highcharts.addEvent(axis.object, "afterSetExtremes", handleAfterSetExtremes)
    Highcharts.addEvent(axis.object.chart, "redraw", handleAfterSetExtremes)

    handleAfterSetExtremes()
    positionInputs()

    // Cleanup: remove ALL three event listeners
    return () => {
      if (axis?.object?.chart && Highcharts) {
        try {
          Highcharts.removeEvent(axis.object.chart, "redraw", positionInputs)
          Highcharts.removeEvent(axis.object, "afterSetExtremes", handleAfterSetExtremes)
          Highcharts.removeEvent(axis.object.chart, "redraw", handleAfterSetExtremes)
        } catch (_e) {
          // Ignore errors during cleanup if objects have been destroyed
        }
      }
    }
  }, [axis, Highcharts, updateChartOptions]) // Added missing dependencies

  const updateRange = () => {
    if (!axis) return
    const extremes = [+minNode.current.value || null, +maxNode.current.value || null]
    axis.setExtremes(...extremes)
  }
  const handleKeyPress = (e) => {
    if (e.key === "Enter") {
      updateRange()
      e.target.select()
    } else {
      if (e.target.name === "min") {
        minNode.current.value = e.target.value
      }
      if (e.target.name === "max") {
        maxNode.current.value = e.target.value
      }
    }
  }

  if (!axis) return null

  const axisExtremes = axis.getExtremes()
  const minColor = !axisExtremes.userMin
    ? "444"
    : axisExtremes.dataMin < axisExtremes.min
      ? "#B1191E"
      : "#5F67CE"
  const maxColor = !axisExtremes.userMax
    ? "444"
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
          placeholder={`${axisExtremes.dataMin || ""}`}
          value={Math.round(min) || ""}
          inputRef={minNode}
          onChange={(e) => updateChartOptions({ extremes: [e.target.value, max] })}
          onKeyPress={handleKeyPress}
          onBlur={updateRange}
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
          placeholder={`${axisExtremes.dataMax || ""}`}
          value={Math.round(max) || ""}
          inputRef={maxNode}
          onChange={(e) => updateChartOptions({ extremes: [min, e.target.value] })}
          onKeyPress={handleKeyPress}
          onBlur={updateRange}
          style={{ ...CLASSES.maxInput, color: maxColor }}
          inputProps={{ style: { textAlign: "center" } }}
        />
      </LightTooltip>
    </div>
  )
}
export default XRangePickers

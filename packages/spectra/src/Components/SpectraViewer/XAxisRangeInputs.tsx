import Input from "@mui/material/Input"
import Tooltip from "@mui/material/Tooltip"
import Highcharts from "highcharts"
import type React from "react"
import { useCallback, useEffect, useRef, useState } from "react"
import { createPortal } from "react-dom"
import { useAxis } from "react-jsx-highcharts"
import { useSpectraStore } from "../../store/spectraStore"

interface XAxisRangeInputsProps {
  enabled?: boolean
  containerId?: string
}

/**
 * Editable x-axis range inputs that replace/overlay the min/max axis labels.
 * Must be rendered as a child of the XAxis component to access the axis instance.
 * Uses a React Portal to render inputs in the correct position.
 */
export const XAxisRangeInputs: React.FC<XAxisRangeInputsProps> = ({
  enabled = true,
  containerId = "spectra-viewer-container",
}) => {
  const axis = useAxis()
  const extremes = useSpectraStore((state) => state.chartOptions.extremes)
  const updateChartOptions = useSpectraStore((state) => state.updateChartOptions)
  const [container, setContainer] = useState<HTMLElement | null>(null)

  // Local state for input values (controlled inputs with validation on blur)
  const [minInput, setMinInput] = useState<string>("")
  const [maxInput, setMaxInput] = useState<string>("")

  // Store extremes in a ref so updateInputsFromAxis doesn't need to depend on it
  const extremesRef = useRef(extremes)

  // Update input values based on current axis state
  // Stable callback (doesn't change when extremes change) to avoid event listener churn
  const updateInputsFromAxis = useCallback(() => {
    if (!axis?.object) return

    const axisMin = axis.object.min
    const axisMax = axis.object.max
    if (axisMin === undefined || axisMax === undefined) return

    const currentExtremes = extremesRef.current

    if (currentExtremes === null) {
      // No extremes set, show actual axis values as integers
      setMinInput(Math.round(axisMin).toString())
      setMaxInput(Math.round(axisMax).toString())
    } else {
      // Extremes are set, show them (or axis values if individual extremes are null)
      const minValue = currentExtremes[0] !== null ? currentExtremes[0] : axisMin
      const maxValue = currentExtremes[1] !== null ? currentExtremes[1] : axisMax
      setMinInput(Math.round(minValue).toString())
      setMaxInput(Math.round(maxValue).toString())
    }
  }, [axis])

  // Update ref and inputs whenever extremes change
  useEffect(() => {
    extremesRef.current = extremes
    updateInputsFromAxis()
  }, [extremes, updateInputsFromAxis])

  // Find the container element
  useEffect(() => {
    const element = document.getElementById(containerId)
    setContainer(element)
  }, [containerId])

  // Position inputs and hide overlapping axis labels
  const [leftPad, setLeftPad] = useState(-5)
  const [rightPad, setRightPad] = useState(0)

  const positionInputsAndHideLabels = useCallback(() => {
    if (!axis?.object) return

    const chart = axis.object.chart
    if (!chart) return

    // Calculate padding based on Y-axis title margins
    let left = -5
    // biome-ignore lint/suspicious/noExplicitAny: Highcharts internal axis properties not typed
    const yAx1 = chart.get("yAx1") as any
    if (yAx1?.axisTitleMargin) {
      left += Number(yAx1.axisTitleMargin)
    }

    let right = 0
    // biome-ignore lint/suspicious/noExplicitAny: Highcharts internal axis properties not typed
    const yAx2 = chart.get("yAx2") as any
    if (yAx2?.axisTitleMargin) {
      right += Number(yAx2.axisTitleMargin)
    }

    setLeftPad(left)
    setRightPad(right)

    // Hide axis labels that overlap with our inputs
    // biome-ignore lint/suspicious/noExplicitAny: Highcharts internal axis properties not typed
    const axisObj = axis.object as any
    if (axisObj.labelGroup?.element?.childNodes) {
      const { min: exMin, max: exMax } = axis.object.getExtremes()
      const tickInterval = axisObj.tickInterval

      axisObj.labelGroup.element.childNodes.forEach((node: Element) => {
        const htmlNode = node as HTMLElement
        // Show all labels first
        htmlNode.style.display = "block"

        // Hide labels that are close to min or max extremes
        const labelValue = Number(htmlNode.textContent)
        const distToMin = Math.abs(labelValue - exMin)
        const distToMax = Math.abs(labelValue - exMax)
        const minDist = Math.min(distToMin, distToMax)

        if (minDist < 0.43 * tickInterval) {
          htmlNode.style.display = "none"
        }
      })
    }
  }, [axis])

  // Listen to chart redraw and update inputs directly
  useEffect(() => {
    if (!axis?.object) return

    const chart = axis.object.chart
    if (!chart) return

    const handleUpdate = () => {
      updateInputsFromAxis()
      positionInputsAndHideLabels()
    }

    // Listen to both redraw and afterSetExtremes to catch all axis changes
    const removeRedraw = Highcharts.addEvent(chart, "redraw", handleUpdate)
    const removeAfterSetExtremes = Highcharts.addEvent(
      axis.object,
      "afterSetExtremes",
      handleUpdate
    )

    // Trigger initial update
    handleUpdate()

    return () => {
      removeRedraw()
      removeAfterSetExtremes()
    }
  }, [axis, updateInputsFromAxis, positionInputsAndHideLabels])

  // Apply extremes from store to chart
  useEffect(() => {
    if (!axis?.object || !enabled) return

    const chart = axis.object.chart

    if (extremes === null) {
      // Reset zoom - Highcharts expects undefined (not null) for full reset
      axis.object.setExtremes(undefined, undefined, true, false)
    } else {
      // Apply extremes (allow individual min/max)
      // undefined means "use default/autoscale", null in our store means "unset"
      const minValue = extremes[0] === null ? undefined : extremes[0]
      const maxValue = extremes[1] === null ? undefined : extremes[1]

      axis.object.setExtremes(
        minValue,
        maxValue,
        true, // redraw
        false // no animation (prevents jarring updates)
      )

      // Only show reset zoom button if at least one extreme is actually set to a number
      // (not just [null, null] which means no zoom)
      const hasActualExtremes = extremes[0] !== null || extremes[1] !== null
      if (hasActualExtremes && chart && typeof chart.showResetZoom === "function") {
        chart.showResetZoom()
      }
    }
  }, [axis, extremes, enabled])

  // Handle min input change
  const handleMinChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    setMinInput(e.target.value)
  }, [])

  // Handle min input blur (validation and commit)
  const handleMinBlur = useCallback(() => {
    if (!minInput.trim()) {
      // Empty input = clear min (keep max)
      updateChartOptions({
        extremes: extremes ? [null, extremes[1]] : null,
      })
      return
    }

    const numValue = parseFloat(minInput)
    if (Number.isNaN(numValue)) {
      // Invalid - reset to current value from extremes or axis
      if (!axis?.object || axis.object.min === undefined) return
      const minValue = extremes?.[0] ?? axis.object.min
      setMinInput(Math.round(minValue).toString())
      return
    }

    // Valid - round to integer and update store
    const roundedValue = Math.round(numValue)
    const currentMax = extremes?.[1] ?? null
    updateChartOptions({ extremes: [roundedValue, currentMax] })
  }, [minInput, extremes, updateChartOptions, axis])

  // Handle max input change
  const handleMaxChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    setMaxInput(e.target.value)
  }, [])

  // Handle max input blur (validation and commit)
  const handleMaxBlur = useCallback(() => {
    if (!maxInput.trim()) {
      // Empty input = clear max (keep min)
      updateChartOptions({
        extremes: extremes ? [extremes[0], null] : null,
      })
      return
    }

    const numValue = parseFloat(maxInput)
    if (Number.isNaN(numValue)) {
      // Invalid - reset to current value from extremes or axis
      if (!axis?.object || axis.object.max === undefined) return
      const maxValue = extremes?.[1] ?? axis.object.max
      setMaxInput(Math.round(maxValue).toString())
      return
    }

    // Valid - round to integer and update store
    const roundedValue = Math.round(numValue)
    const currentMin = extremes?.[0] ?? null
    updateChartOptions({ extremes: [currentMin, roundedValue] })
  }, [maxInput, extremes, updateChartOptions, axis])

  // Handle Enter key to commit input
  const handleKeyDown = useCallback((e: React.KeyboardEvent) => {
    if (e.key === "Enter") {
      ;(e.target as HTMLInputElement).blur()
    }
  }, [])

  if (!enabled || !container) return null

  const inputsElement = (
    <div
      style={{
        height: 0,
        position: "relative",
        bottom: 38,
        pointerEvents: "none",
      }}
    >
      <Tooltip
        title="Type to change min, clear to autoscale"
        placement="top-end"
        slotProps={{
          transition: { timeout: { enter: 150, exit: 400 } },
          tooltip: {
            sx: {
              backgroundColor: "white",
              color: "rgba(0, 0, 0, 0.87)",
              boxShadow: 1,
              fontSize: 12,
              margin: "0 13px 7px",
            },
          },
        }}
      >
        <Input
          value={minInput}
          onChange={handleMinChange}
          onBlur={handleMinBlur}
          onKeyDown={handleKeyDown}
          type="text"
          sx={{
            position: "absolute",
            left: `${leftPad}px`,
            width: 30,
            fontWeight: "bold",
            fontSize: "0.75rem",
            color: "#444",
            pointerEvents: "auto",
            "& input": {
              textAlign: "center",
            },
          }}
        />
      </Tooltip>

      <Tooltip
        title="Type to change max, clear to autoscale"
        placement="top-start"
        slotProps={{
          transition: { timeout: { enter: 150, exit: 400 } },
          tooltip: {
            sx: {
              backgroundColor: "white",
              color: "rgba(0, 0, 0, 0.87)",
              boxShadow: 1,
              fontSize: 12,
              margin: "0 13px 7px",
            },
          },
        }}
      >
        <Input
          type="text"
          value={maxInput}
          onChange={handleMaxChange}
          onBlur={handleMaxBlur}
          onKeyDown={handleKeyDown}
          sx={{
            position: "absolute",
            right: `${rightPad}px`,
            width: 30,
            fontWeight: "bold",
            fontSize: "0.75rem",
            color: "#444",
            pointerEvents: "auto",
            "& input": {
              textAlign: "center",
            },
          }}
        />
      </Tooltip>
    </div>
  )

  return createPortal(inputsElement, container)
}

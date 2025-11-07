import Input from "@mui/material/Input"
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

  // Listen to chart redraw and update inputs directly
  useEffect(() => {
    if (!axis?.object) return

    const chart = axis.object.chart
    if (!chart) return

    // Listen to both redraw and afterSetExtremes to catch all axis changes
    const removeRedraw = Highcharts.addEvent(chart, "redraw", updateInputsFromAxis)
    const removeAfterSetExtremes = Highcharts.addEvent(
      axis.object,
      "afterSetExtremes",
      updateInputsFromAxis
    )

    // Trigger initial update
    updateInputsFromAxis()

    return () => {
      removeRedraw()
      removeAfterSetExtremes()
    }
  }, [axis, updateInputsFromAxis])

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
        position: "absolute",
        bottom: 2,
        left: 10,
        right: 10,
        display: "flex",
        alignItems: "center",
        justifyContent: "flex-start",
        gap: 8,
        zIndex: 100,
        pointerEvents: "none", // Allow chart interactions through the div
      }}
    >
      <Input
        value={minInput}
        onChange={handleMinChange}
        onBlur={handleMinBlur}
        onKeyDown={handleKeyDown}
        size="small"
        type="text"
        inputProps={{
          step: 1,
          min: 300,
          max: 1000,
        }}
        sx={{
          width: 80,
          pointerEvents: "auto", // Re-enable for input
          "& input": {
            fontSize: "0.75rem",
            padding: "4px 8px",
            textAlign: "center",
          },
          "& .MuiOutlinedInput-root": {
            backgroundColor: "rgba(255, 255, 255, 0.9)",
          },
        }}
      />

      <Input
        type="text"
        value={maxInput}
        onChange={handleMaxChange}
        onBlur={handleMaxBlur}
        onKeyDown={handleKeyDown}
        size="small"
        inputProps={{
          step: 1,
          min: 300,
          max: 1000,
        }}
        sx={{
          width: 80,
          pointerEvents: "auto",
          "& input": {
            fontSize: "0.75rem",
            padding: "4px 8px",
            textAlign: "center",
          },
          "& .MuiOutlinedInput-root": {
            backgroundColor: "rgba(255, 255, 255, 0.9)",
          },
        }}
      />
    </div>
  )

  return createPortal(inputsElement, container)
}

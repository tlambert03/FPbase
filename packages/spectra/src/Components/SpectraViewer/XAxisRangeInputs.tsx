import Input from "@mui/material/Input"
import Tooltip from "@mui/material/Tooltip"
import Highcharts from "highcharts"
import type React from "react"
import { useCallback, useEffect, useMemo, useRef, useState } from "react"
import { createPortal } from "react-dom"
import { useAxis } from "react-jsx-highcharts"
import { useSpectraStore } from "../../store/spectraStore"

interface XAxisRangeInputsProps {
  enabled?: boolean
  containerId?: string
  extremes?: [number | null, number | null] | null
}

/**
 * Editable x-axis range inputs that replace/overlay the min/max axis labels.
 * Must be rendered as a child of the XAxis component to access the axis instance.
 * Uses a React Portal to render inputs in the correct position.
 */
export const XAxisRangeInputs: React.FC<XAxisRangeInputsProps> = ({
  enabled = true,
  containerId = "spectra-viewer-container",
  extremes: extremesProp,
}) => {
  const axis = useAxis()
  const storeExtremes = useSpectraStore((state) => state.chartOptions.extremes)
  const updateChartOptions = useSpectraStore((state) => state.updateChartOptions)

  // Use provided extremes prop if available (for SimpleSpectraViewer with provideState)
  // Otherwise fall back to store extremes (for main app)
  const extremes = extremesProp !== undefined ? extremesProp : storeExtremes
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
  const prevPaddingRef = useRef({ left: -5, right: 0 })

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

    // Only update state if values changed
    if (left !== prevPaddingRef.current.left) {
      setLeftPad(left)
      prevPaddingRef.current.left = left
    }
    if (right !== prevPaddingRef.current.right) {
      setRightPad(right)
      prevPaddingRef.current.right = right
    }

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
    if (!axis?.object) return

    const chart = axis.object.chart

    if (extremes === null) {
      // Reset zoom - Highcharts expects undefined (not null) for full reset
      if (enabled) {
        axis.object.setExtremes(undefined, undefined, true, false)
      }
      // Always hide reset zoom button when extremes are null, even if component is disabled
      // This ensures the button disappears when "Remove All Spectra" is clicked
      // biome-ignore lint/suspicious/noExplicitAny: Highcharts resetZoomButton not typed
      if (chart && (chart as any).resetZoomButton) {
        // biome-ignore lint/suspicious/noExplicitAny: Highcharts resetZoomButton not typed
        ;(chart as any).resetZoomButton = (chart as any).resetZoomButton.destroy()
      }
    } else if (enabled) {
      // Only apply extremes if component is enabled
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

  // Shared blur handler logic for min/max inputs
  const handleRangeBlur = useCallback(
    (
      isMin: boolean,
      inputValue: string,
      setInputValue: (value: string) => void,
      axisValue: number | undefined
    ) => {
      if (!inputValue.trim()) {
        // Empty input = clear this extreme (keep the other)
        const newExtremes: [number | null, number | null] | null = extremes
          ? isMin
            ? [null, extremes[1]]
            : [extremes[0], null]
          : null
        updateChartOptions({ extremes: newExtremes })
        return
      }

      const numValue = parseFloat(inputValue)
      if (Number.isNaN(numValue)) {
        // Invalid - reset to current value from extremes or axis
        if (axisValue === undefined) return
        const fallbackValue = isMin ? (extremes?.[0] ?? axisValue) : (extremes?.[1] ?? axisValue)
        setInputValue(Math.round(fallbackValue).toString())
        return
      }

      // Valid - round to integer and update store
      const roundedValue = Math.round(numValue)
      const newExtremes: [number | null, number | null] = isMin
        ? [roundedValue, extremes?.[1] ?? null]
        : [extremes?.[0] ?? null, roundedValue]
      updateChartOptions({ extremes: newExtremes })
    },
    [extremes, updateChartOptions]
  )

  // Handle min input change
  const handleMinChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    setMinInput(e.target.value)
  }, [])

  // Handle min input blur (validation and commit)
  const handleMinBlur = useCallback(() => {
    handleRangeBlur(true, minInput, setMinInput, axis?.object?.min)
  }, [handleRangeBlur, minInput, axis])

  // Handle max input change
  const handleMaxChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    setMaxInput(e.target.value)
  }, [])

  // Handle max input blur (validation and commit)
  const handleMaxBlur = useCallback(() => {
    handleRangeBlur(false, maxInput, setMaxInput, axis?.object?.max)
  }, [handleRangeBlur, maxInput, axis])

  // Handle Enter key to commit input
  const handleKeyDown = useCallback((e: React.KeyboardEvent) => {
    if (e.key === "Enter") {
      ;(e.target as HTMLInputElement).blur()
    }
  }, [])

  // Calculate input colors based on clipping state
  // Gray = autoscale, Red = clipping data, Blue = user zoom (no clipping)
  // IMPORTANT: Must be before early return to satisfy Rules of Hooks
  const minColor = useMemo(() => {
    if (!extremes || extremes[0] === null) return "#444" // Autoscale
    if (!axis?.object) return "#444"
    // biome-ignore lint/suspicious/noExplicitAny: dataMin/dataMax not in Highcharts type definitions
    const axisObj = axis.object as any
    const isClipping =
      axisObj.dataMin !== undefined &&
      axis.object.min !== undefined &&
      axisObj.dataMin < axis.object.min
    return isClipping ? "#B1191E" : "#5F67CE"
  }, [extremes, axis])

  const maxColor = useMemo(() => {
    if (!extremes || extremes[1] === null) return "#444" // Autoscale
    if (!axis?.object) return "#444"
    // biome-ignore lint/suspicious/noExplicitAny: dataMin/dataMax not in Highcharts type definitions
    const axisObj = axis.object as any
    const isClipping =
      axisObj.dataMax !== undefined &&
      axis.object.max !== undefined &&
      axisObj.dataMax > axis.object.max
    return isClipping ? "#B1191E" : "#5F67CE"
  }, [extremes, axis])

  if (!enabled || !container) return null

  // Shared tooltip configuration
  const tooltipSlotProps = {
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
  }

  // Shared input styles (position-specific styles applied inline)
  const baseInputSx = {
    position: "absolute" as const,
    width: 30,
    fontWeight: "bold",
    fontSize: "0.75rem",
    pointerEvents: "auto" as const,
    "& input": {
      textAlign: "center" as const,
    },
  }

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
        slotProps={tooltipSlotProps}
      >
        <Input
          name="min"
          value={minInput}
          onChange={handleMinChange}
          onBlur={handleMinBlur}
          onKeyDown={handleKeyDown}
          type="text"
          sx={{ ...baseInputSx, left: `${leftPad}px`, color: minColor }}
        />
      </Tooltip>

      <Tooltip
        title="Type to change max, clear to autoscale"
        placement="top-start"
        slotProps={tooltipSlotProps}
      >
        <Input
          type="text"
          name="max"
          value={maxInput}
          onChange={handleMaxChange}
          onBlur={handleMaxBlur}
          onKeyDown={handleKeyDown}
          sx={{ ...baseInputSx, right: `${rightPad}px`, color: maxColor }}
        />
      </Tooltip>
    </div>
  )

  return createPortal(inputsElement, container)
}

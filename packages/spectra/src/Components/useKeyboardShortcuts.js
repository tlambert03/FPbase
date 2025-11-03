import { useCallback, useEffect } from "react"
import { useSpectraStore } from "../store/spectraStore"

// Palette order for cycling (must match keys in palettes.js)
const PALETTE_ORDER = ["wavelength", "rainbow", "tol_contrast", "tol_vibrant", "okabe_ito"]

/**
 * Get the next palette in the cycle
 * @param {string} currentPalette - Current palette name
 * @returns {string} Next palette name
 */
const getNextPalette = (currentPalette) => {
  const currentIndex = PALETTE_ORDER.indexOf(currentPalette)
  return PALETTE_ORDER[(currentIndex + 1) % PALETTE_ORDER.length]
}

const useKeyboardShortcuts = () => {
  const chartOptions = useSpectraStore((state) => state.chartOptions)
  const updateChartOptions = useSpectraStore((state) => state.updateChartOptions)

  const toggleY = useCallback(
    () => updateChartOptions({ showY: !chartOptions.showY }),
    [chartOptions.showY, updateChartOptions]
  )
  const toggleX = useCallback(
    () => updateChartOptions({ showX: !chartOptions.showX }),
    [chartOptions.showX, updateChartOptions]
  )
  const toggleGrid = useCallback(
    () => updateChartOptions({ showGrid: !chartOptions.showGrid }),
    [chartOptions.showGrid, updateChartOptions]
  )
  const toggleScaleEC = useCallback(
    () => updateChartOptions({ scaleEC: !chartOptions.scaleEC }),
    [chartOptions.scaleEC, updateChartOptions]
  )
  const toggleScaleQY = useCallback(
    () => updateChartOptions({ scaleQY: !chartOptions.scaleQY }),
    [chartOptions.scaleQY, updateChartOptions]
  )
  const toggleShareTooltip = useCallback(
    () => updateChartOptions({ shareTooltip: !chartOptions.shareTooltip }),
    [chartOptions.shareTooltip, updateChartOptions]
  )
  const toggleAreaFill = useCallback(
    () => updateChartOptions({ areaFill: !chartOptions.areaFill }),
    [chartOptions.areaFill, updateChartOptions]
  )

  const cyclePalette = useCallback(() => {
    const nextPalette = getNextPalette(chartOptions.palette || "wavelength")
    updateChartOptions({ palette: nextPalette })
  }, [chartOptions.palette, updateChartOptions])

  useEffect(() => {
    const handleKeyDown = (event) => {
      // don't do anything if we're on an input
      if (document.activeElement && document.activeElement.tagName.toUpperCase() === "INPUT") {
        return
      }
      switch (event.code) {
        case "KeyX":
          event.preventDefault()
          toggleX()
          break
        case "KeyY":
          event.preventDefault()
          toggleY()
          break
        case "KeyG":
          event.preventDefault()
          toggleGrid()
          break
        case "KeyE":
          event.preventDefault()
          toggleScaleEC()
          break
        case "KeyQ":
          event.preventDefault()
          toggleScaleQY()
          break
        case "KeyA":
        case "KeyF":
          event.preventDefault()
          toggleAreaFill()
          break
        case "KeyT":
          event.preventDefault()
          toggleShareTooltip()
          break
        case "KeyC":
          event.preventDefault()
          cyclePalette()
          break
        default:
          break
      }
    }

    document.addEventListener("keydown", handleKeyDown)
    return () => {
      document.removeEventListener("keydown", handleKeyDown)
    }
  }, [
    toggleAreaFill,
    toggleGrid,
    toggleScaleEC,
    toggleScaleQY,
    toggleShareTooltip,
    toggleX,
    toggleY,
    cyclePalette,
  ])
}

export default useKeyboardShortcuts

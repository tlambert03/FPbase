import { memo, useMemo } from "react"
import { Series } from "react-jsx-highcharts"
import { useSpectrum } from "../../hooks/useSpectraQueries"
import PALETTES from "../../palettes"

const OD = (num) => (num <= 0 ? 10 : -Math.log10(num))

const hex2rgba = (hex, alpha = 1) => {
  if (!hex) return `rgba(0,0,0,${alpha})` // Default to black if no color provided
  const [r, g, b] = hex.match(/\w\w/g).map((x) => parseInt(x, 16))
  return `rgba(${r},${g},${b},${alpha})`
}

const CROSS_HATCH = {
  pattern: {
    path: {
      d: ["M 5,5 L 10,10", "M 5,5 L 0,10", "M 5,5 L 10,0", "M 5,5 L 0,0"],
    },
    width: 10,
    height: 10,
    color: "#ddd",
    opacity: 0.4,
  },
}

const VERT_LINES = {
  pattern: {
    path: {
      d: ["M 2,10 L 2,0"],
    },
    width: 3,
    height: 10,
    opacity: 0.2,
  },
}

// ErrorBoundary removed - was non-functional (all logic commented out)
// If error handling is needed in the future, implement a proper error boundary

/**
 * Hook that applies all data transformations to spectrum data
 * Ensures all transformations are computed from the original spectrum.data
 * with proper memoization dependencies, making transformations invertible
 */
const useTransformedSpectrumData = ({
  spectrum,
  ownerInfo,
  exNorm,
  scaleEC,
  scaleQY,
  inverted,
  logScale,
}) => {
  // Determine if we need to fetch EX spectrum for normalization
  const needsExNorm = (spectrum.subtype === "EM" || spectrum.subtype === "O") && exNorm

  // Find the EX spectrum ID if needed
  const exSpectrumId = useMemo(() => {
    if (!needsExNorm || !spectrum.owner?.slug || !(spectrum.owner.slug in ownerInfo)) return null

    const ownerSpectra = ownerInfo[spectrum.owner.slug].spectra
    if (!ownerSpectra) return null

    const exSpectrum =
      ownerSpectra.find((i) => i.subtype === "EX") || ownerSpectra.find((i) => i.subtype === "AB")

    return exSpectrum ? exSpectrum.id : null
  }, [needsExNorm, spectrum.owner?.slug, ownerInfo])

  // Fetch the EX spectrum data if needed
  const { data: exSpectrumData } = useSpectrum(exSpectrumId)

  // Determine if transformations should be applied
  const willScaleEC = Boolean(
    (spectrum.subtype === "EX" || spectrum.subtype === "AB") && scaleEC && spectrum.owner?.extCoeff
  )
  const willScaleQY = Boolean(
    (spectrum.subtype === "EM" || spectrum.subtype === "O") && scaleQY && spectrum.owner?.qy
  )

  // Apply ALL transformations in a single memoized pipeline
  // This ensures transformations are always computed from the original spectrum.data,
  // making them properly invertible when toggled
  const transformedData = useMemo(() => {
    // Start fresh from original data
    let data = [...spectrum.data]

    // 1. ExNorm transformation (if applicable)
    if (needsExNorm && exSpectrumData) {
      // Find the scalar value at the exNorm wavelength
      let scalar = 0
      const exEfficiency = exSpectrumData.data.find(([x]) => x === exNorm)
      if (exEfficiency) {
        ;[, scalar] = exEfficiency
      }
      data = data.map(([a, b]) => [a, b * scalar])
    }

    // 2. Extinction Coefficient scaling (for EX/AB spectra)
    if (willScaleEC) {
      data = data.map(([a, b]) => [a, b * spectrum.owner.extCoeff])
    }

    // 3. Quantum Yield scaling (for EM/O spectra)
    if (willScaleQY) {
      data = data.map(([a, b]) => [a, b * spectrum.owner.qy])
    }

    // 4. Inversion transformation
    if (inverted) {
      data = data.map(([a, b]) => [a, 1 - b])
    }

    // 5. Log scale transformation
    if (logScale) {
      data = data.map(([a, b]) => [a, OD(b)])
    }

    return data
  }, [
    spectrum.data,
    spectrum.owner?.extCoeff,
    spectrum.owner?.qy,
    needsExNorm,
    exSpectrumData,
    exNorm,
    willScaleEC,
    willScaleQY,
    inverted,
    logScale,
  ])

  return transformedData
}

const SpectrumSeries = memo(function SpectrumSeries({
  inverted,
  logScale,
  scaleEC,
  scaleQY,
  spectrum,
  areaFill,
  exNorm,
  palette,
  ownerIndex,
  ownerInfo,
  visible = true,
}) {
  // All transformations are now handled in the hook with proper memoization
  // Note: Hook must be called before any conditional returns (Rules of Hooks)
  const serie = useTransformedSpectrumData({
    spectrum,
    ownerInfo,
    exNorm,
    scaleEC,
    scaleQY,
    inverted,
    logScale,
  })

  if (!spectrum) return null

  // Determine if scaling was applied (for display purposes)
  const willScaleEC = Boolean(
    (spectrum.subtype === "EX" || spectrum.subtype === "AB") && scaleEC && spectrum.owner?.extCoeff
  )
  const willScaleQY = Boolean(
    (spectrum.subtype === "EM" || spectrum.subtype === "O") && scaleQY && spectrum.owner?.qy
  )

  let name = `${spectrum.owner.name}`
  if (["EX", "EM", "2P", "AB"].includes(spectrum.subtype)) {
    name += ` ${spectrum.subtype}`
  }
  let dashStyle = "Solid"
  if (["EX", "AB", "2P"].includes(spectrum.subtype)) {
    dashStyle = "ShortDash"
  }
  let myColor = spectrum.color
  if (palette !== "wavelength" && palette in PALETTES) {
    const { hexlist } = PALETTES[palette]
    myColor = hexlist[ownerIndex % hexlist.length]
  }
  let color = hex2rgba(myColor, 0.9)
  let fillColor = hex2rgba(myColor, 0.5)
  let lineWidth = areaFill ? 0.5 : 1.8
  let type = areaFill ? "areaspline" : "spline"
  if (spectrum.category === "C") {
    fillColor = CROSS_HATCH
  }
  if (spectrum.category === "L") {
    fillColor = { ...VERT_LINES }
    lineWidth = areaFill ? 1 : 1.8
  }
  if (["BS", "LP"].includes(spectrum.subtype)) {
    lineWidth = 1.8
    type = "spline"
    color = "#999"
  }

  return (
    <Series
      type={type}
      subtype={spectrum.subtype}
      scaleEC={willScaleEC}
      scaleQY={willScaleQY}
      name={name}
      visible={visible}
      color={color}
      fillColor={fillColor}
      dashStyle={dashStyle}
      lineWidth={lineWidth}
      className={`cat-${spectrum.category} subtype-${spectrum.subtype}`}
      data={serie}
      threshold={logScale ? 10 : 0}
    />
  )
})

export default SpectrumSeries

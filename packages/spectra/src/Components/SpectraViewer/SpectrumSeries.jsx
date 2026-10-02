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

const isEx = (spectrum) => spectrum.subtype === "EX" || spectrum.subtype === "AB"
const isEm = (spectrum) => spectrum.subtype === "EM" || spectrum.subtype === "O"

/** ID of the owner's EX (or AB) spectrum, needed to normalize an emission spectrum to exNorm */
export function exNormSpectrumId(spectrum, ownerInfo, exNorm) {
  if (!(isEm(spectrum) && exNorm)) return null
  const ownerSpectra = ownerInfo?.[spectrum.owner?.slug]?.spectra
  if (!ownerSpectra) return null
  const ex =
    ownerSpectra.find((i) => i.subtype === "EX") || ownerSpectra.find((i) => i.subtype === "AB")
  return ex ? ex.id : null
}

/**
 * Highcharts series options for one spectrum.
 *
 * All transformations are applied to the original `spectrum.data`, so toggling them is
 * invertible.
 */
export function spectrumSeriesOptions(
  spectrum,
  { inverted, logScale, scaleEC, scaleQY, areaFill, exNorm, palette },
  { exSpectrum, ownerIndex, visible, yAxis }
) {
  const willScaleEC = Boolean(isEx(spectrum) && scaleEC && spectrum.owner?.extCoeff)
  const willScaleQY = Boolean(isEm(spectrum) && scaleQY && spectrum.owner?.qy)

  let data = [...spectrum.data]
  if (isEm(spectrum) && exNorm && exSpectrum) {
    const exEfficiency = exSpectrum.data.find(([x]) => x === exNorm)
    const scalar = exEfficiency ? exEfficiency[1] : 0
    data = data.map(([a, b]) => [a, b * scalar])
  }
  if (willScaleEC) data = data.map(([a, b]) => [a, b * spectrum.owner.extCoeff])
  if (willScaleQY) data = data.map(([a, b]) => [a, b * spectrum.owner.qy])
  if (inverted) data = data.map(([a, b]) => [a, 1 - b])
  if (logScale) data = data.map(([a, b]) => [a, OD(b)])

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

  return {
    id: String(spectrum.id),
    yAxis,
    type,
    subtype: spectrum.subtype,
    scaleEC: willScaleEC,
    scaleQY: willScaleQY,
    name,
    visible,
    color,
    fillColor,
    dashStyle,
    lineWidth,
    className: `cat-${spectrum.category} subtype-${spectrum.subtype}`,
    data,
    threshold: logScale ? 10 : 0,
  }
}

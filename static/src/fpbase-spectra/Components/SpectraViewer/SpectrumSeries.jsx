import React, { memo } from "react"
import { Series } from "react-jsx-highcharts"

const OD = num => (num <= 0 ? 10 : -Math.log10(num))

const hex2rgba = (hex, alpha = 1) => {
  const [r, g, b] = hex.match(/\w\w/g).map(x => parseInt(x, 16));
  return `rgba(${r},${g},${b},${alpha})`;
};

const CROSS_HATCH = {
  pattern: {
    path: {
      d: ["M 5,5 L 10,10", "M 5,5 L 0,10", "M 5,5 L 10,0", "M 5,5 L 0,0"]
    },
    width: 10,
    height: 10,
    color: "#ddd",
    opacity: 0.4
  }
}

const VERT_LINES = {
  pattern: {
    path: {
      d: "M 2,10 L 2,0"
    },
    width: 3,
    height: 10,
    opacity: 0.2
  }
}

const SpectrumSeries = memo(function SpectrumSeries({
  inverted,
  logScale,
  scaleEC,
  scaleQY,
  spectrum,
  areaFill
}) {
  if (!spectrum) return
  const willScaleEC = Boolean(
    spectrum.subtype === "EX" && scaleEC && spectrum.owner.extCoeff
  )
  const willScaleQY = Boolean(
    spectrum.subtype === "EM" && scaleQY && spectrum.owner.qy
  )

  let serie = [...spectrum.data]
  if (willScaleEC) {
    serie = serie.map(([a, b]) => [a, b * spectrum.owner.extCoeff])
  }
  if (spectrum.subtype === "EM") {
    if (scaleQY && spectrum.owner.qy) {
      serie = serie.map(([a, b]) => [a, b * spectrum.owner.qy])
    }
  }
  if (inverted) {
    serie = serie.map(([a, b]) => [a, 1 - b])
  }
  if (logScale) {
    serie = [...serie].map(([a, b]) => [a, OD(b)])
  }

  let name = `${spectrum.owner.name}`
  if (["EX", "EM", "A_2P", "2P"].includes(spectrum.subtype)) {
    name += ` ${spectrum.subtype.replace("A_", "")}`
  }
  let color = hex2rgba(spectrum.color, 0.9)
  let fillColor = hex2rgba(spectrum.color, 0.5)
  let lineWidth = 0.5
  let type = areaFill ? "areaspline" : 'spline'
  if (spectrum.category === "C") {
    fillColor = CROSS_HATCH
  }
  if (spectrum.category === "L") {
    fillColor = {...VERT_LINES}
    lineWidth = 1
  }
  if (["BS", "LP"].includes(spectrum.subtype)) {
    lineWidth = 1.5
    type="spline"
    color="#999"
  }
  return (
    <Series
      type={type}
      subtype={spectrum.subtype}
      scaleEC={willScaleEC}
      scaleQY={willScaleQY}
      name={name}
      color={color}
      fillColor={fillColor}
      lineWidth={lineWidth}
      className={`cat-${spectrum.category} subtype-${spectrum.subtype}`}
      data={serie}
      threshold={logScale ? 10 : 0}
    />
  )
})

export default SpectrumSeries

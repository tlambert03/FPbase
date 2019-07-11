import qs from "qs"

function stateToUrl(activeSpectra, chartOptions, exNorm) {
  const qstrings = []
  if (activeSpectra) qstrings.push(`s=${activeSpectra.join(",")}`)
  if (chartOptions) {
    const opts = { ...chartOptions }
    delete opts.__typename
    const [xMin, xMax] = opts.extremes
    if (xMin) opts.xMin = xMin
    if (xMax) opts.xMax = xMax
    delete opts.extremes
    Object.keys(opts).forEach(key => {
      if (typeof opts[key] === "boolean") {
        opts[key] = Number(opts[key])
      }
    })
    qstrings.push(qs.stringify(opts))
  }
  if (exNorm[0] && exNorm[1]) {
    qstrings.push(qs.stringify({ normWave: exNorm[0], normID: exNorm[1] }))
  }
  const { origin, pathname } = window.location
  return qstrings.length ? `${origin + pathname}?${qstrings.join("&")}` : ""
}

export default stateToUrl

import qs from "qs"
import { defaults } from "./defaults"
import PALETTES from "./palettes"

// import { GET_SPECTRUM } from "./queries"
// import client from "./client"

// const DEFAULT_EXPIRY = 0 * 60 * 60 // n hours

// const ID = () => {
//   // Math.random should be unique because of its seeding algorithm.
//   // Convert it to base 36 (numbers + letters), and grab the first 9 characters
//   // after the decimal.
//   return `_${Math.random()
//     .toString(36)
//     .substr(2, 9)}`
// }

// const emptyFormSelector = () => ({ id: ID(), value: null })

const getStorageWithExpire = (cacheKey, expiry = 12 * 60 * 60) => {
  const cached = localStorage.getItem(cacheKey)
  const whenCached = localStorage.getItem(`${cacheKey}:ts`)
  if (cached !== null && whenCached !== null) {
    const age = (Date.now() - whenCached) / 1000
    if (age < expiry) {
      return JSON.parse(cached)
    }
    // We need to clean up this old key
    localStorage.removeItem(cacheKey)
    localStorage.removeItem(`${cacheKey}:ts`)
  }
  return null
}

const setStorageWithTimeStamp = (cacheKey, value) => {
  localStorage.setItem(cacheKey, JSON.stringify(value))
  localStorage.setItem(`${cacheKey}:ts`, Date.now())
}

// const getCachedSpectrum = async (id, options) => {
//   let expiry = DEFAULT_EXPIRY // 10 min default
//   if (typeof options === "number") {
//     expiry = options
//   } else if (typeof options === "object") {
//     // I hope you didn't set it to 0 seconds
//     expiry = options.seconds || expiry
//   }

//   const cacheKey = `spectrum_${id}`
//   const cached = getStorageWithExpire(cacheKey, expiry)
//   if (cached !== null) return Promise.resolve(JSON.parse(cached))

//   const { loading, error, data } = await client.query({
//     query: GET_SPECTRUM,
//     variables: { id }
//   })
//   if (!loading && !error) {
//     const spec = data.spectrum
//     localStorage.setItem(cacheKey, JSON.stringify(spec))
//     localStorage.setItem(`${cacheKey}:ts`, Date.now())
//     return Promise.resolve(spec)
//   }
//   return Promise.reject()
// }

// async function fetchSpectrum(id) {
//   const data = await getCachedSpectrum(id)
//   if (data) {
//     return {
//       ...data,
//       name: `${data.owner.name}${
//         ["P", "D"].includes(data.category)
//           ? ` ${data.subtype.replace(/^A_/g, "")}`
//           : ""
//       }`,
//       inverted: false,
//       visible: true,
//       ecNormed: false,
//       qyNormed: false
//     }
//   }
//   const msg = `Could not find spectrum with ID: ${id}`
//   throw new Error(msg)
// }

// async function fetchSpectraList(ids) {
//   const results = await Promise.all(
//     ids.map(id => fetchSpectrum(id)).map(p => p.catch(e => e))
//   )
//   const valid = results.filter(result => !(result instanceof Error))
//   return valid
// }

const debounce = (fn, time) => {
  let timeout

  return function (...args) {
    const functionCall = () => fn.apply(this, args)

    clearTimeout(timeout)
    timeout = setTimeout(functionCall, time)
  }
}

function reshapeSpectraInfo(arr) {
  if (!arr) return {}
  return arr.reduce(
    (prev, cur) => {
      if (!Object.hasOwn(prev.ownerInfo, cur.owner.slug)) {
        // eslint-disable-next-line no-param-reassign
        prev.ownerInfo[cur.owner.slug] = {
          category: cur.category.toUpperCase(),
          label: cur.owner.name,
          spectra: [],
          value: cur.owner.slug,
          url: cur.owner.url,
        }
      }
      prev.ownerInfo[cur.owner.slug].spectra.push({
        id: String(cur.id),
        subtype: cur.subtype.toUpperCase(),
        active: true,
      })
      // eslint-disable-next-line no-param-reassign
      prev.spectraInfo[cur.id] = {
        subtype: cur.subtype.toUpperCase(),
        owner: cur.owner.slug,
        label: cur.owner.name,
        category: cur.category.toUpperCase(),
      }
      return prev
    },
    { ownerInfo: {}, spectraInfo: {} }
  )
}

function decoder(str, _, charset) {
  const strWithoutPlus = str.replace(/\+/g, " ")
  if (charset === "iso-8859-1") {
    // unescape never throws, no try...catch needed:
    return strWithoutPlus.replace(/%[0-9a-f]{2}/gi, unescape)
  }

  if (/^(\d+|\d*\.\d+)$/.test(str)) {
    return str
    // return parseFloat(str)
  }

  const keywords = {
    true: true,
    false: false,
    null: null,
    undefined,
  }
  if (str in keywords) {
    return keywords[str]
  }

  // utf-8
  try {
    return decodeURIComponent(strWithoutPlus)
  } catch (_e) {
    return strWithoutPlus
  }
}

function isTouchDevice() {
  try {
    const prefixes = " -webkit- -moz- -o- -ms- ".split(" ")

    const mq = (query) => window.matchMedia(query).matches

    if (
      "ontouchstart" in window ||
      (typeof window.DocumentTouch !== "undefined" && document instanceof window.DocumentTouch)
    ) {
      return true
    }

    return mq(["(", prefixes.join("touch-enabled),("), "heartz", ")"].join(""))
  } catch (_e) {
    // console.error("(Touch detect failed)", e)
    return false
  }
}

function trapz(arr) {
  // approximate area under curve as series of trapezoids
  // arr = [[wave, data], ...]
  let sum = 0
  for (let i = 1; i < arr.length; i++) {
    sum += 0.5 * (arr[i][1] + arr[i - 1][1]) * (arr[i][0] - arr[i - 1][0])
  }
  return sum
}

function spectraProduct(ar1, ar2) {
  // calculate product of two spectra.values
  // these assume monotonic increase w/ step = 1
  const output = []
  const left = Math.max(ar1[0][0], ar2[0][0]) // the min wavelength shared by both arrays
  const right = Math.min(ar1[ar1.length - 1][0], ar2[ar2.length - 1][0]) // the max wavelength shared by both arrays
  if (left >= right) return []

  const offsetA1 = left - ar1[0][0]
  const offsetA2 = left - ar2[0][0]
  for (let i = 0; i < right - left; i++) {
    output.push([left + i, ar1[offsetA1 + i][1] * ar2[offsetA2 + i][1]])
  }
  return output
}

/**
 * Check if an ID is valid
 * @param {string|number} id - The ID to validate
 * @returns {boolean} Whether the ID is valid
 */
const isValidId = (id) => {
  if (!id) return false
  if (!Number.isNaN(parseFloat(id))) return true
  if (typeof id === "string") {
    if (id.startsWith("$cl") || id.startsWith("$cf")) {
      return true
    }
    return id.split("_").every((i) => isValidId(i))
  }
  return false
}

/**
 * Filter spectra IDs to only valid ones
 * @param {string[]} spectra - Array of spectrum IDs
 * @returns {string[]} Array of valid spectrum IDs
 */
export const validSpectraIds = (spectra) => spectra.filter((id) => isValidId(id))

/**
 * Parse URL parameters to extract chart options and active spectra
 * @param {Object} data - Initial data object (defaults to defaults from defaults.js)
 * @returns {Object} Parsed data with activeSpectra, chartOptions, exNorm, etc.
 */
export function parseURL(data) {
  data = data || defaults
  const url = qs.parse(window.location.search.replace(/^\?/, ""), {
    decoder,
  })
  if (Object.keys(url).length === 0 && url.constructor === Object) return data

  const booleanOptions = Object.keys(defaults.chartOptions).filter(
    (key) => typeof defaults.chartOptions[key] === "boolean"
  )

  const extremes = [null, null]
  const exNorm = [null, null]
  Object.keys(url).forEach((key) => {
    if (booleanOptions.includes(key)) {
      data.chartOptions[key] = Boolean(+url[key])
    }
    if (key === "palette" && url[key] in PALETTES) {
      data.chartOptions.palette = url[key]
    }
    if (key === "xMin") extremes[0] = +url[key]
    if (key === "xMax") extremes[1] = +url[key]
    if (key === "normWave") exNorm[0] = url[key]
    if (key === "normID") exNorm[1] = url[key]
    if (["s", "activeSpectra"].includes(key)) {
      let active = url[key]
      if (!Array.isArray(active)) active = active.split(",")
      data.activeSpectra = validSpectraIds(active)
    }
  })
  if (extremes.some((i) => i)) data.chartOptions.extremes = extremes
  if (exNorm.some((i) => i)) data.exNorm = exNorm
  return data
}

export {
  debounce,
  reshapeSpectraInfo,
  decoder,
  getStorageWithExpire,
  setStorageWithTimeStamp,
  isTouchDevice,
  trapz,
  spectraProduct,
}

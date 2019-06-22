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

  return function() {
    const functionCall = () => fn.apply(this, arguments)

    clearTimeout(timeout)
    timeout = setTimeout(functionCall, time)
  }
}

function reshapeSpectraInfo(arr) {
  if (!arr) return {}
  return arr.reduce(
    (prev, cur) => {
      if (!Object.prototype.hasOwnProperty.call(prev.owners, cur.owner.slug)) {
        // eslint-disable-next-line no-param-reassign
        prev.owners[cur.owner.slug] = {
          category: cur.category,
          label: cur.owner.name,
          spectra: [],
          value: cur.owner.slug,
          url: cur.owner.url
        }
      }
      prev.owners[cur.owner.slug].spectra.push({
        id: cur.id,
        subtype: cur.subtype,
        active: true
      })
      // eslint-disable-next-line no-param-reassign
      prev.spectraInfo[cur.id] = {
        subtype: cur.subtype,
        owner: cur.owner.slug,
        label: cur.owner.name,
        category: cur.category
      }
      return prev
    },
    { owners: {}, spectraInfo: {} }
  )
}

function decoder(str, decoder, charset) {
  const strWithoutPlus = str.replace(/\+/g, " ")
  if (charset === "iso-8859-1") {
    // unescape never throws, no try...catch needed:
    return strWithoutPlus.replace(/%[0-9a-f]{2}/gi, unescape)
  }

  if (/^(\d+|\d*\.\d+)$/.test(str)) {
    return str
    //return parseFloat(str)
  }

  const keywords = {
    true: true,
    false: false,
    null: null,
    undefined
  }
  if (str in keywords) {
    return keywords[str]
  }

  // utf-8
  try {
    return decodeURIComponent(strWithoutPlus)
  } catch (e) {
    return strWithoutPlus
  }
}

export {
  debounce,
  reshapeSpectraInfo,
  decoder,
  getStorageWithExpire,
  setStorageWithTimeStamp
}

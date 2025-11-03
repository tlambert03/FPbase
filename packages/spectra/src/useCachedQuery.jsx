import { useEffect, useState } from "react"
import "unfetch/polyfill/index"

// Local helper functions (only used by this legacy hook)
const getStorageWithExpire = (cacheKey, expiry = 12 * 60 * 60) => {
  const cached = localStorage.getItem(cacheKey)
  const whenCached = localStorage.getItem(`${cacheKey}:ts`)
  if (cached !== null && whenCached !== null) {
    const age = (Date.now() - whenCached) / 1000
    if (age < expiry) {
      return JSON.parse(cached)
    }
    localStorage.removeItem(cacheKey)
    localStorage.removeItem(`${cacheKey}:ts`)
  }
  return null
}

const setStorageWithTimeStamp = (cacheKey, value) => {
  localStorage.setItem(cacheKey, JSON.stringify(value))
  localStorage.setItem(`${cacheKey}:ts`, Date.now())
}

const useCachedFetch = (url, cacheKey, maxAge) => {
  // fetch data no more than once every maxAge
  const [stash, setStash] = useState(getStorageWithExpire(cacheKey, maxAge))
  useEffect(() => {
    if (!stash) {
      fetch(url)
        .then((i) => i.json())
        .then(({ data }) => {
          setStash(data)
          setStorageWithTimeStamp(cacheKey, data)
        })
    }
  }, [cacheKey, stash, url]) // eslint-disable-line
  return stash
}

export { useCachedFetch }

import { useEffect, useState } from "react"

/**
 * Retrieves data from localStorage with an expiration check.
 *
 * @param cacheKey - The localStorage key to retrieve
 * @param expiry - Cache lifetime in seconds (default: 12 hours)
 * @returns Cached data or null if expired/missing
 */
const getStorageWithExpire = <T>(cacheKey: string, expiry = 12 * 60 * 60): T | null => {
  const cached = localStorage.getItem(cacheKey)
  const whenCached = localStorage.getItem(`${cacheKey}:ts`)
  if (cached !== null && whenCached !== null) {
    const age = (Date.now() - Number(whenCached)) / 1000
    if (age < expiry) {
      return JSON.parse(cached) as T
    }
    localStorage.removeItem(cacheKey)
    localStorage.removeItem(`${cacheKey}:ts`)
  }
  return null
}

/**
 * Stores data in localStorage with a timestamp.
 *
 * @param cacheKey - The localStorage key to store under
 * @param value - The data to cache
 */
const setStorageWithTimeStamp = <T>(cacheKey: string, value: T): void => {
  localStorage.setItem(cacheKey, JSON.stringify(value))
  localStorage.setItem(`${cacheKey}:ts`, String(Date.now()))
}

/**
 * Custom hook that fetches data with localStorage caching and expiration.
 * Legacy hook used for simple data fetching with client-side caching.
 *
 * @param url - The URL to fetch data from
 * @param cacheKey - localStorage key for caching
 * @param maxAge - Cache lifetime in seconds
 * @returns Cached or fetched data, or null while loading
 */
export const useCachedFetch = <T>(url: string, cacheKey: string, maxAge: number): T | null => {
  // fetch data no more than once every maxAge
  const [stash, setStash] = useState<T | null>(getStorageWithExpire<T>(cacheKey, maxAge))

  useEffect(() => {
    if (!stash) {
      fetch(url)
        .then((i) => i.json())
        .then(({ data }: { data: T }) => {
          setStash(data)
          setStorageWithTimeStamp(cacheKey, data)
        })
    }
  }, [cacheKey, stash, url]) // eslint-disable-line

  return stash
}

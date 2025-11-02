import { useEffect, useState } from "react"
import { getStorageWithExpire, setStorageWithTimeStamp } from "./util"
import "unfetch/polyfill/index"

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

import { useState, useEffect } from "react"
import { useApolloClient } from "react-apollo-hooks"
import { getStorageWithExpire, setStorageWithTimeStamp } from "./util"

export const useCachedQuery = (query, cacheKey, maxAge) => {
  // fetch data no more than once every maxAge
  const [stash, setStash] = useState(getStorageWithExpire(cacheKey, maxAge))
  const client = useApolloClient()
  useEffect(() => {
    async function fetchData() {
      const { data } = await client.query({ query })
      if (data) {
        setStash(data)
        setStorageWithTimeStamp(cacheKey, data)
      }
    }
    if (!stash) fetchData()
  }, []) // eslint-disable-line
  return stash
}

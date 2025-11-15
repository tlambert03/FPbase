import { useQuery } from "@tanstack/react-query"
import { useEffect } from "react"
import { useMetadataStore } from "../store/metadataStore"
import type { SpectraListResponse } from "../types"
import { reshapeSpectraInfo } from "../utils/spectraUtils"

/**
 * Fetch and cache spectra metadata from REST endpoint
 * Uses ETag-based caching for efficient updates
 */
export function useSpectraMetadata() {
  const { setMetadata } = useMetadataStore()

  const { data, isLoading, error } = useQuery({
    queryKey: ["spectraList"],
    queryFn: async () => {
      const response = await fetch("/api/spectra-list/")
      if (!response.ok) {
        throw new Error(`Failed to fetch spectra: ${response.statusText}`)
      }
      const json = await response.json()
      const data = json.data as SpectraListResponse
      return data.spectra
    },
    staleTime: 5 * 60 * 1000, // 5 minutes - time after which data is considered stale
    gcTime: 60 * 60 * 1000, // 1 hour - time after which unused cache is garbage collected
    refetchOnWindowFocus: false, // Don't refetch on focus
    refetchOnReconnect: false, // Don't refetch on reconnect
  })

  // Update metadata store when data changes
  useEffect(() => {
    if (data) {
      const { ownerInfo, spectraInfo } = reshapeSpectraInfo(data)
      setMetadata(ownerInfo, spectraInfo)
    }
  }, [data, setMetadata])

  return { data, isLoading, error }
}

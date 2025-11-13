import { useQuery } from "@tanstack/react-query"
import { useEffect } from "react"
import { fetchGraphQL } from "../api/client"
import type { SpectraListResponse } from "../api/queries"
import { SPECTRA_LIST } from "../api/queries"
import { useMetadataStore } from "../store/metadataStore"
import { reshapeSpectraInfo } from "../utils/spectraUtils"

/**
 * Fetch and cache spectra metadata using GraphQL
 * Migrated from REST endpoint to use GraphQL with server-side caching
 */
export function useSpectraMetadata() {
  const { setMetadata } = useMetadataStore()

  const { data, isLoading, error } = useQuery({
    queryKey: ["spectraList"],
    queryFn: async () =>
      (await fetchGraphQL<SpectraListResponse>(SPECTRA_LIST, { method: "GET" })).spectra,
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

import { useQuery } from "@tanstack/react-query"
import { useEffect } from "react"
import { fetchGraphQL } from "../api/client"
import { SPECTRA_LIST } from "../api/queries"
import { useMetadataStore } from "../store/metadataStore"
import { reshapeSpectraInfo } from "../utils/spectraUtils"

interface SpectraSlug {
  id: string
  category: string
  subtype: string
  owner: {
    id: string
    name: string
    slug: string
    url: string | null
  }
}

interface SpectraListResponse {
  spectra: SpectraSlug[]
}

/**
 * Fetch and cache spectra metadata using GraphQL
 * Migrated from REST endpoint to use GraphQL with server-side caching
 */
export function useSpectraMetadata() {
  const { setMetadata } = useMetadataStore()

  const { data, isLoading, error } = useQuery({
    queryKey: ["spectra-list"],
    queryFn: async () => {
      const response = await fetchGraphQL<SpectraListResponse>(SPECTRA_LIST, undefined, {
        operationName: "_FPB_SpectraList",
      })
      return response.spectra
    },
    staleTime: 10 * 60 * 1000, // 10 minutes
    gcTime: 30 * 60 * 1000, // 30 minutes (formerly cacheTime)
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

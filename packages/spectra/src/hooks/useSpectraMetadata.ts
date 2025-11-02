import { useQuery } from "@tanstack/react-query"
import { useEffect } from "react"
import { fetchAPI } from "../api/client"
import { useMetadataStore } from "../store/metadataStore"
import { reshapeSpectraInfo } from "../utils/spectraUtils"

interface SpectraSlug {
  id: string
  category: string
  subtype: string
  owner: {
    name: string
    slug: string
    url: string
  }
}

interface SpectraSlugsResponse {
  data: {
    spectra: SpectraSlug[]
  }
}

/**
 * Fetch and cache spectra metadata
 * Replaces useCachedFetch with TanStack Query
 */
export function useSpectraMetadata() {
  const { setMetadata } = useMetadataStore()

  const { data, isLoading, error } = useQuery({
    queryKey: ["spectra-slugs"],
    queryFn: async () => {
      const response = await fetchAPI<SpectraSlugsResponse>("/api/proteins/spectraslugs/")
      return response.data.spectra
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

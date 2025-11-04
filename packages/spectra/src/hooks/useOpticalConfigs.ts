import { useQuery } from "@tanstack/react-query"
import { fetchGraphQL } from "../api/client"
import { OPTICAL_CONFIG_LIST } from "../api/queries"

interface Microscope {
  id: number
  name: string
}

interface OpticalConfig {
  id: number
  name: string
  comments: string
  microscope: Microscope
}

interface OpticalConfigsResponse {
  opticalConfigs: OpticalConfig[]
}

/**
 * Fetch and cache optical configurations list
 * Replaces useCachedFetch for optical configs with TanStack Query
 */
export function useOpticalConfigs() {
  return useQuery({
    queryKey: ["opticalConfigs"],
    queryFn: async (): Promise<OpticalConfig[]> => {
      const response = await fetchGraphQL<OpticalConfigsResponse>(OPTICAL_CONFIG_LIST)
      return response.opticalConfigs
    },
    staleTime: 10 * 60 * 1000, // 10 minutes (matches old cache duration)
    gcTime: 30 * 60 * 1000, // 30 minutes garbage collection
  })
}

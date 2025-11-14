import { useQuery } from "@tanstack/react-query"

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

interface OpticalConfigsListResponse {
  opticalConfigs: OpticalConfig[]
}

/**
 * Fetch and cache optical configurations list from REST endpoint
 * Uses ETag-based caching for efficient updates
 */
export function useOpticalConfigs() {
  return useQuery({
    queryKey: ["opticalConfigs"],
    queryFn: async (): Promise<OpticalConfig[]> => {
      const response = await fetch("/api/optical-configs-list/")
      if (!response.ok) {
        throw new Error(`Failed to fetch optical configs: ${response.statusText}`)
      }
      const json = await response.json()
      const data = json.data as OpticalConfigsListResponse
      return data.opticalConfigs
    },
    staleTime: 10 * 60 * 1000, // 10 minutes (matches cache duration)
    gcTime: 30 * 60 * 1000, // 30 minutes garbage collection
    refetchOnWindowFocus: false, // Don't refetch on focus
    refetchOnReconnect: false, // Don't refetch on reconnect
  })
}

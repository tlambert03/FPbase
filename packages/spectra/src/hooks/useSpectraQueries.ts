import { useQueries, useQuery } from "@tanstack/react-query"
import { fetchGraphQL } from "../api/client"
import { GET_OPTICAL_CONFIG, GET_SPECTRUM } from "../api/queries"
import type { OpticalConfig, Spectrum } from "../types"

/**
 * Normalize spectrum subtype from API
 * Graphene-django prefixes "A_" to enum values that aren't valid GraphQL identifiers
 * (e.g., "2P" starts with a digit, so becomes "A_2P"). We normalize to the canonical "2P".
 */
function normalizeSpectrum<T extends Spectrum | null>(spectrum: T): T {
  if (!spectrum) return spectrum
  return {
    ...spectrum,
    subtype: spectrum.subtype === "A_2P" ? "2P" : spectrum.subtype,
  } as T
}

/**
 * Fetch a single spectrum by ID
 */
export function useSpectrum(id: string | null) {
  return useQuery({
    queryKey: ["spectrum", id],
    queryFn: async () => {
      if (!id) return null
      const data = await fetchGraphQL<{ spectrum: Spectrum }>(GET_SPECTRUM, {
        id: Number.parseInt(id, 10),
      })
      return data.spectrum
    },
    select: normalizeSpectrum,
    enabled: !!id,
    staleTime: 5 * 60 * 1000, // 5 minutes
  })
}

/**
 * Batch fetch multiple spectra by IDs using individual queries
 * This approach provides better caching and handles empty arrays correctly
 */
export function useSpectraBatch(ids: string[]) {
  return useQueries({
    queries: ids.map((id) => ({
      queryKey: ["spectrum", id],
      queryFn: async () => {
        const data = await fetchGraphQL<{ spectrum: Spectrum }>(GET_SPECTRUM, {
          id: Number.parseInt(id, 10),
        })
        return data.spectrum
      },
      select: normalizeSpectrum,
      staleTime: 5 * 60 * 1000, // 5 minutes
    })),
    combine: (results) => {
      // Extract all non-null spectrum data
      const data = results.map((result) => result.data).filter((s): s is Spectrum => s !== null)

      // Return a stable object shape matching useQuery's return type
      // This ensures proper re-renders when the data changes
      return {
        data,
        isPending: results.some((result) => result.isPending),
        isError: results.some((result) => result.isError),
        error: results.find((result) => result.error)?.error ?? null,
      }
    },
  })
}

/**
 * Fetch optical configuration by ID
 */
export function useOpticalConfig(id: number | null) {
  return useQuery({
    queryKey: ["optical-config", id],
    queryFn: async () => {
      if (!id) return null
      const data = await fetchGraphQL<{ opticalConfig: OpticalConfig }>(GET_OPTICAL_CONFIG, { id })
      return data.opticalConfig
    },
    enabled: !!id,
    staleTime: 10 * 60 * 1000, // 10 minutes
  })
}

/**
 * TanStack Query configuration
 * Use this when setting up QueryClientProvider
 */
export const queryClientConfig = {
  defaultOptions: {
    queries: {
      retry: 1,
      refetchOnWindowFocus: false,
      staleTime: 5 * 60 * 1000, // 5 minutes default
    },
  },
}

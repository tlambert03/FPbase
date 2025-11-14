import { useQuery } from "@tanstack/react-query"
import { fetchAPI, fetchGraphQL } from "../api/client"
import { batchSpectraQuery, GET_OPTICAL_CONFIG, GET_SPECTRUM, SPECTRA_LIST } from "../api/queries"
import type { OpticalConfig, OpticalConfigInfo, SpectraListResponse, Spectrum } from "../types"

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
 * Normalize subtype for any object with a subtype field
 * Used for list responses that don't include full spectrum data
 */
function normalizeSubtype<T extends { subtype: string }>(item: T): T {
  return {
    ...item,
    subtype: item.subtype === "A_2P" ? "2P" : item.subtype,
  }
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
 * Batch fetch multiple spectra by IDs
 * More efficient than individual queries
 */
export function useSpectraBatch(ids: string[]) {
  return useQuery({
    queryKey: ["spectra-batch", ids.sort().join(",")],
    queryFn: async () => {
      if (ids.length === 0) return []

      const query = batchSpectraQuery(ids)
      const data = await fetchGraphQL<Record<string, Spectrum>>(query, {})

      // Convert object to array
      return Object.values(data)
    },
    select: (data) => data.map(normalizeSpectrum),
    enabled: ids.length > 0,
    staleTime: 5 * 60 * 1000, // 5 minutes
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
 * Fetch list of all spectra (for search/autocomplete)
 */
export function useSpectraList() {
  return useQuery({
    queryKey: ["spectra-list"],
    queryFn: async () => {
      const data = await fetchGraphQL<SpectraListResponse>(SPECTRA_LIST, {})
      return data.spectra
    },
    select: (data) => data.map(normalizeSubtype),
    staleTime: 10 * 60 * 1000, // 10 minutes
  })
}

/**
 * Fetch optical configs info from REST API
 * Uses REST instead of GraphQL for better caching
 */
export function useOpticalConfigsInfo() {
  return useQuery({
    queryKey: ["optical-configs-info"],
    queryFn: async () => {
      const data = await fetchAPI<{ opticalConfigs: OpticalConfigInfo[] }>("/api/proteins/ocinfo/")
      return data.opticalConfigs
    },
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

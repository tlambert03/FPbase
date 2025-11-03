import { useMemo } from "react"
import COLORS from "../colors"
import { useSpectraStore } from "../store/spectraStore"
import type { Spectrum, SpectrumSubtype } from "../types"
import { computeOverlap, trapz } from "../utils/spectraUtils"
import { useSpectraBatch } from "./useSpectraQueries"

const rangexy = (start: number, end: number): number[] =>
  Array.from({ length: end - start }, (_v, k) => k + start)

/**
 * Generate custom laser spectrum (e.g., $cl1_488)
 * Format: $cl<id>_<wavelength>
 */
function generateCustomLaser(id: string): Spectrum {
  const parts = id.split("_")
  const customId = parts[0] ?? id
  const waveStr = parts[1] ?? "500"
  const wave = Number.parseInt(waveStr, 10)

  const data: [number, number][] = [
    [wave - 1, 0],
    [wave, 1],
    [wave + 1, 0],
  ]

  return {
    id: customId,
    customId: id,
    subtype: "PD",
    owner: {
      id: id,
      slug: id,
      name: `${wave} laser`,
    },
    category: "L",
    data,
    area: trapz(data),
    color: wave in COLORS ? COLORS[wave] : "#999999",
  }
}

/**
 * Generate custom filter spectrum (e.g., $cf1_bp_500_50_90)
 * Format: $cf<id>_<type>_<center>_<width>_<transmission>
 */
function generateCustomFilter(id: string): Spectrum {
  const parts = id.split("_")
  const customId = parts[0] ?? id
  const subtypeStr = parts[1] ?? "BP"
  const centerStr = parts[2] ?? "500"
  const widthStr = parts[3] ?? "50"
  const transStr = parts[4]

  const subtype = subtypeStr.toUpperCase() as "BP" | "LP" | "SP"
  const center = Number.parseInt(centerStr, 10)
  const width = Number.parseInt(widthStr, 10)
  const trans = transStr ? Number.parseInt(transStr, 10) / 100 : 0.9

  const data: [number, number][] = []
  let name = "Custom "

  switch (subtype) {
    case "BP": {
      // Band pass filter
      const min = Math.round(center - width / 2)
      const max = Math.round(center + width / 2)
      data.push([min - 1, 0])
      for (const x of rangexy(min, max + 1)) {
        data.push([x, trans])
      }
      data.push([max + 1, 0])
      name += `${center}/${width} bp`
      break
    }
    case "LP":
      // Long pass filter
      for (const x of rangexy(300, center)) {
        data.push([x, 0])
      }
      for (const x of rangexy(center + 1, 1000)) {
        data.push([x, trans])
      }
      name += `${center}lp`
      break
    case "SP":
      // Short pass filter
      for (const x of rangexy(300, center)) {
        data.push([x, trans])
      }
      for (const x of rangexy(center + 1, 1000)) {
        data.push([x, 0])
      }
      name += `${center}sp`
      break
  }

  return {
    id: customId,
    customId: id,
    subtype: subtype as SpectrumSubtype,
    owner: {
      id,
      slug: id,
      name,
    },
    category: "F",
    data,
    area: trapz(data),
    color: center in COLORS ? COLORS[center] : "#999999",
  }
}

/**
 * Pure function to compute visible spectra from inputs
 * Separated for testability and clarity
 */
function computeVisibleSpectra(
  activeSpectra: string[],
  activeOverlaps: string[],
  apiSpectra: Spectrum[] | undefined,
  overlapCache: Record<string, Spectrum>,
  hiddenSpectra: string[],
  isProvidedMode: boolean,
  setOverlapCache: (id: string, spectrum: Spectrum) => void
): Spectrum[] {
  // Separate custom and real spectrum IDs
  const customIds = activeSpectra.filter((id) => id?.startsWith("$c"))

  // Generate custom spectra
  const customSpectra: Spectrum[] = customIds
    .map((id) => {
      if (id.startsWith("$cf")) return generateCustomFilter(id)
      if (id.startsWith("$cl")) return generateCustomLaser(id)
      return null
    })
    .filter((s): s is Spectrum => s !== null)

  // Add area calculation to API spectra
  const apiSpectraWithArea: Spectrum[] = (apiSpectra || []).map((s) => ({
    ...s,
    area: s.area ?? trapz(s.data),
  }))

  // Combine available spectra for overlap computation
  const allAvailableSpectra = [...apiSpectraWithArea, ...customSpectra]

  // Compute overlaps from activeOverlaps IDs
  const overlapSpectra: Spectrum[] = activeOverlaps
    .map((id) => {
      // Check cache first
      if (overlapCache[id]) return overlapCache[id]

      // Parse overlap ID to get constituent spectrum IDs (e.g., "123_456")
      const spectrumIds = id.split("_")
      const constituents = spectrumIds
        .map((sid) => allAvailableSpectra.find((s) => (s.customId || s.id) === sid))
        .filter((s): s is Spectrum => !!s)
        .filter((s): s is Spectrum & { owner: NonNullable<Spectrum["owner"]> } => s.owner !== null)

      // Only compute if we have all constituent spectra with valid owners
      if (constituents.length === spectrumIds.length) {
        const overlap = computeOverlap(...constituents)
        setOverlapCache(id, overlap)
        return overlap
      }

      return null
    })
    .filter((s): s is Spectrum => !!s)

  // Combine all spectra
  const allSpectra = [...allAvailableSpectra, ...overlapSpectra]

  // Filter out hidden spectra (only when using store data, not provided IDs)
  return isProvidedMode
    ? allSpectra
    : allSpectra.filter((s) => !hiddenSpectra.includes(s.customId || s.id))
}

/**
 * Hook to fetch and manage spectra data
 * Handles both real spectra from API and custom generated spectra
 *
 * @param providedIds - Optional array of spectrum IDs to use (overrides store)
 * @param providedOverlaps - Optional array of overlap IDs to use (overrides store)
 * @returns Array of spectrum objects with data
 */
export function useSpectraData(providedIds?: string[], providedOverlaps?: string[]): Spectrum[] {
  // Get active spectra from store if not provided
  const storeActiveSpectra = useSpectraStore((state) => state.activeSpectra)
  const storeActiveOverlaps = useSpectraStore((state) => state.activeOverlaps)
  const hiddenSpectra = useSpectraStore((state) => state.hiddenSpectra)
  const overlapCache = useSpectraStore((state) => state.overlapCache)
  const setOverlapCache = useSpectraStore((state) => state.setOverlapCache)

  const activeSpectra = providedIds ?? storeActiveSpectra
  const activeOverlaps = providedOverlaps ?? storeActiveOverlaps

  // Separate real IDs for batch fetching
  const realIds = useMemo(
    () => activeSpectra.filter((id) => id && !id.startsWith("$c")),
    [activeSpectra]
  )

  // Fetch real spectra using TanStack Query
  const { data: apiSpectra } = useSpectraBatch(realIds)

  // Compute final spectrum list (pure derivation, no side effects)
  return useMemo(
    () =>
      computeVisibleSpectra(
        activeSpectra,
        activeOverlaps,
        apiSpectra,
        overlapCache,
        hiddenSpectra,
        providedIds !== undefined, // Fix: check undefined, not truthiness
        setOverlapCache
      ),
    [
      activeSpectra,
      activeOverlaps,
      apiSpectra,
      overlapCache,
      hiddenSpectra,
      providedIds,
      setOverlapCache,
    ]
  )
}

export default useSpectraData

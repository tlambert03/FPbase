import { useEffect, useState } from "react"
import COLORS from "../colors"
import { useOverlapCache } from "../store/overlapsStore"
import { useSpectraStore } from "../store/spectraStore"
import type { Spectrum, SpectrumSubtype } from "../types"
import { useSpectraBatch } from "./useSpectraQueries"

const rangexy = (start: number, end: number): number[] =>
  Array.from({ length: end - start }, (_v, k) => k + start)

// Verify TypeScript version is loading
console.log("✅ NEW TypeScript useSpectraData.ts loaded!")
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
    subtype: "pd",
    owner: {
      id: id,
      slug: id,
      name: `${wave} laser`,
    },
    category: "L",
    data,
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
    color: center in COLORS ? COLORS[center] : "#999999",
  }
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
  const [currentData, setCurrentData] = useState<Spectrum[]>([])

  // Get active spectra from store if not provided
  const storeActiveSpectra = useSpectraStore((state) => state.activeSpectra)
  const storeActiveOverlaps = useSpectraStore((state) => state.activeOverlaps)
  const overlapCache = useOverlapCache()

  const activeSpectra = providedIds ?? storeActiveSpectra
  const activeOverlaps = providedOverlaps ?? storeActiveOverlaps

  // Separate custom and real spectrum IDs
  const customIds = activeSpectra.filter((id) => id?.startsWith("$c"))
  const realIds = activeSpectra.filter((id) => id && !id.startsWith("$c"))

  // Fetch real spectra using TanStack Query
  const { data: apiSpectra } = useSpectraBatch(realIds)

  useEffect(() => {
    // Generate custom spectra
    const customSpectra: Spectrum[] = customIds
      .map((id) => {
        if (id.startsWith("$cf")) {
          return generateCustomFilter(id)
        }
        if (id.startsWith("$cl")) {
          return generateCustomLaser(id)
        }
        return null
      })
      .filter((s): s is Spectrum => s !== null)

    // Get overlap data from cache
    const overlapSpectra: Spectrum[] = activeOverlaps
      .map((id) => overlapCache[id])
      .filter((s): s is Spectrum => !!s)

    // Combine all spectra
    const allSpectra = [...(apiSpectra || []), ...customSpectra, ...overlapSpectra]

    // Only update if data actually changed
    const currentIds = currentData.map((s) => s.customId || s.id).join(",")
    const newIds = allSpectra.map((s) => s.customId || s.id).join(",")

    if (currentIds !== newIds) {
      setCurrentData(allSpectra)
    }
  }, [activeOverlaps, apiSpectra, customIds, currentData, overlapCache])

  return currentData
}

export default useSpectraData

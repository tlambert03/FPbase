import { create } from "zustand"
import type { Spectrum } from "../types"

/**
 * Overlaps store for caching computed spectrum overlaps
 * Replaces window.OverlapCache
 */

interface OverlapsState {
  cache: Record<string, Spectrum>
  getOverlap: (id: string) => Spectrum | undefined
  setOverlap: (id: string, spectrum: Spectrum) => void
  clearCache: () => void
}

export const useOverlapsStore = create<OverlapsState>((set, get) => ({
  cache: {},

  getOverlap: (id) => get().cache[id],

  setOverlap: (id, spectrum) =>
    set((state) => ({
      cache: { ...state.cache, [id]: spectrum },
    })),

  clearCache: () => set({ cache: {} }),
}))

// Selector hook
export const useOverlapCache = () => useOverlapsStore((state) => state.cache)

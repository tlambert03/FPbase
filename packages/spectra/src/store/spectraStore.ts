import { create } from "zustand"
import { createJSONStorage, persist } from "zustand/middleware"
import { defaults } from "../defaults"
import type { SpectraStore } from "../types"

// Helper function to clear URL params when user modifies state
function clearUrlIfNeeded(get: () => SpectraStore) {
  const state = get()
  if (state._urlInitialized && window.location.search) {
    // User is modifying state - clear URL to show it's no longer valid
    window.history.replaceState({}, "", window.location.pathname)
    // Return true to indicate the URL was cleared
    return true
  }
  return false
}

export const useSpectraStore = create<SpectraStore>()(
  persist(
    (set, get) => ({
      // Initial state
      activeSpectra: defaults.activeSpectra,
      activeOverlaps: defaults.activeOverlaps,
      hiddenSpectra: [],
      excludeSubtypes: defaults.excludeSubtypes,
      exNorm: defaults.exNorm,
      chartOptions: defaults.chartOptions,
      customFilters: {},
      customLasers: {},
      overlapCache: {}, // Computed data, not persisted (see partialize)
      _urlInitialized: false,

      // Spectra management
      setActiveSpectra: (ids) => {
        const wasCleared = clearUrlIfNeeded(get)
        set({ activeSpectra: ids, _urlInitialized: wasCleared ? false : get()._urlInitialized })
      },

      updateActiveSpectra: (add = [], remove = []) => {
        const wasCleared = clearUrlIfNeeded(get)
        set((state) => {
          let newSpectra = [...state.activeSpectra]

          // Remove spectra
          if (remove.length > 0) {
            newSpectra = newSpectra.filter((id) => !remove.includes(id))
          }

          // Add spectra (avoid duplicates)
          if (add.length > 0) {
            const toAdd = add.filter((id) => !newSpectra.includes(id))
            newSpectra = [...newSpectra, ...toAdd]
          }

          return {
            activeSpectra: newSpectra,
            _urlInitialized: wasCleared ? false : state._urlInitialized,
          }
        })
      },

      // Overlaps management
      setActiveOverlaps: (ids) => {
        const wasCleared = clearUrlIfNeeded(get)
        set({ activeOverlaps: ids, _urlInitialized: wasCleared ? false : get()._urlInitialized })
      },

      updateActiveOverlaps: (add = [], remove = []) => {
        const wasCleared = clearUrlIfNeeded(get)
        set((state) => {
          let newOverlaps = [...state.activeOverlaps]

          // Remove overlaps
          if (remove.length > 0) {
            newOverlaps = newOverlaps.filter((id) => !remove.includes(id))
          }

          // Add overlaps (avoid duplicates)
          if (add.length > 0) {
            const toAdd = add.filter((id) => !newOverlaps.includes(id))
            newOverlaps = [...newOverlaps, ...toAdd]
          }

          return {
            activeOverlaps: newOverlaps,
            _urlInitialized: wasCleared ? false : state._urlInitialized,
          }
        })
      },

      // Visibility management
      toggleSpectrumVisibility: (id) => {
        const wasCleared = clearUrlIfNeeded(get)
        set((state) => ({
          hiddenSpectra: state.hiddenSpectra.includes(id)
            ? state.hiddenSpectra.filter((hid) => hid !== id)
            : [...state.hiddenSpectra, id],
          _urlInitialized: wasCleared ? false : state._urlInitialized,
        }))
      },

      setHiddenSpectra: (ids) => {
        const wasCleared = clearUrlIfNeeded(get)
        set({ hiddenSpectra: ids, _urlInitialized: wasCleared ? false : get()._urlInitialized })
      },

      // Subtype management
      setExcludeSubtypes: (subtypes) => set({ excludeSubtypes: subtypes }),

      // Excitation normalization
      setExNorm: (norm) => {
        const wasCleared = clearUrlIfNeeded(get)
        set({ exNorm: norm, _urlInitialized: wasCleared ? false : get()._urlInitialized })
      },

      // Chart options
      updateChartOptions: (options) => {
        const wasCleared = clearUrlIfNeeded(get)
        set((state) => ({
          chartOptions: { ...state.chartOptions, ...options },
          _urlInitialized: wasCleared ? false : state._urlInitialized,
        }))
      },

      // Custom spectra
      addCustomFilter: (filter) =>
        set((state) => ({
          customFilters: {
            ...state.customFilters,
            [filter.id]: filter,
          },
        })),

      removeCustomFilter: (id) =>
        set((state) => {
          const { [id]: _removed, ...rest } = state.customFilters
          return { customFilters: rest }
        }),

      addCustomLaser: (laser) =>
        set((state) => ({
          customLasers: {
            ...state.customLasers,
            [laser.id]: laser,
          },
        })),

      removeCustomLaser: (id) =>
        set((state) => {
          const { [id]: _removed, ...rest } = state.customLasers
          return { customLasers: rest }
        }),

      // Overlap cache management
      setOverlapCache: (id, spectrum) =>
        set((state) => ({
          overlapCache: { ...state.overlapCache, [id]: spectrum },
        })),

      clearOverlapCache: () => set({ overlapCache: {} }),

      // URL initialization tracking
      setUrlInitialized: (value) => set({ _urlInitialized: value }),
    }),
    {
      name: "fpbase-spectra-storage",
      version: 1, // Increment this to invalidate old stored state
      storage: createJSONStorage(() => sessionStorage),
      // Only persist certain fields (_urlInitialized is intentionally excluded)
      partialize: (state) => ({
        activeSpectra: state.activeSpectra,
        hiddenSpectra: state.hiddenSpectra,
        excludeSubtypes: state.excludeSubtypes,
        chartOptions: state.chartOptions,
        customFilters: state.customFilters,
        customLasers: state.customLasers,
      }),
    }
  )
)

// Selector hooks for optimized re-renders
export const useActiveSpectra = () => useSpectraStore((state) => state.activeSpectra)
export const useChartOptions = () => useSpectraStore((state) => state.chartOptions)
export const useExcludeSubtypes = () => useSpectraStore((state) => state.excludeSubtypes)

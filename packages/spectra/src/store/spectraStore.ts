import { create } from "zustand"
import { createJSONStorage, persist } from "zustand/middleware"
import { defaults } from "../defaults"
import type { SpectraStore } from "../types"

export const useSpectraStore = create<SpectraStore>()(
  persist(
    (set) => ({
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
      setActiveSpectra: (ids) => set({ activeSpectra: ids }),

      updateActiveSpectra: (add = [], remove = []) =>
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

          // Auto-cleanup: clear exNorm if removing the normed spectrum
          const newExNorm =
            remove.length > 0 && state.exNorm?.[1] && remove.includes(state.exNorm[1])
              ? defaults.exNorm
              : state.exNorm

          return {
            activeSpectra: newSpectra,
            exNorm: newExNorm,
          }
        }),

      // Overlaps management
      setActiveOverlaps: (ids) => set({ activeOverlaps: ids }),

      updateActiveOverlaps: (add = [], remove = []) =>
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

          return { activeOverlaps: newOverlaps }
        }),

      // Visibility management
      toggleSpectrumVisibility: (id) =>
        set((state) => ({
          hiddenSpectra: state.hiddenSpectra.includes(id)
            ? state.hiddenSpectra.filter((hid) => hid !== id)
            : [...state.hiddenSpectra, id],
        })),

      setHiddenSpectra: (ids) => set({ hiddenSpectra: ids }),

      // Subtype management
      setExcludeSubtypes: (subtypes) => set({ excludeSubtypes: subtypes }),

      // Excitation normalization
      setExNorm: (norm) => {
        // Don't clear URL for exNorm changes - it's part of URL state
        set({ exNorm: norm })
      },

      // Chart options
      updateChartOptions: (options) => {
        // Don't clear URL for chart option changes - they're part of URL state
        set((state) => ({
          chartOptions: { ...state.chartOptions, ...options },
        }))
      },

      // Custom spectra management
      addCustomFilter: (id, params) =>
        set((state) => ({
          customFilters: {
            ...state.customFilters,
            [id]: params,
          },
        })),

      updateCustomFilter: (id, params) =>
        set((state) => {
          const existing = state.customFilters[id]
          if (!existing) return state
          return {
            customFilters: {
              ...state.customFilters,
              [id]: { ...existing, ...params },
            },
          }
        }),

      removeCustomFilter: (id) =>
        set((state) => {
          const { [id]: _removed, ...rest } = state.customFilters
          return { customFilters: rest }
        }),

      addCustomLaser: (id, params) =>
        set((state) => ({
          customLasers: {
            ...state.customLasers,
            [id]: params,
          },
        })),

      updateCustomLaser: (id, params) =>
        set((state) => {
          const existing = state.customLasers[id]
          if (!existing) return state
          return {
            customLasers: {
              ...state.customLasers,
              [id]: { ...existing, ...params },
            },
          }
        }),

      removeCustomLaser: (id) =>
        set((state) => {
          const { [id]: _removed, ...rest } = state.customLasers

          // Auto-cleanup: clear exNorm if it references this laser
          const newExNorm = state.exNorm?.[1] === id ? defaults.exNorm : state.exNorm

          return {
            customLasers: rest,
            exNorm: newExNorm,
          }
        }),

      // Overlap cache management
      setOverlapCache: (id, spectrum) =>
        set((state) => ({
          overlapCache: { ...state.overlapCache, [id]: spectrum },
        })),

      clearOverlapCache: () => set({ overlapCache: {} }),

      // Clear all spectra and related state
      clearAllSpectra: () =>
        set((state) => ({
          activeSpectra: [],
          activeOverlaps: [],
          hiddenSpectra: [],
          exNorm: defaults.exNorm,
          chartOptions: {
            ...state.chartOptions,
            extremes: null,
          },
          customFilters: {},
          customLasers: {},
          overlapCache: {},
        })),

      // URL initialization tracking
      setUrlInitialized: (value) => set({ _urlInitialized: value }),

      // Atomically replace state (single render, no mixing with existing state)
      replace: (state) =>
        set({
          activeSpectra: state.activeSpectra ?? defaults.activeSpectra,
          activeOverlaps: state.activeOverlaps ?? defaults.activeOverlaps,
          hiddenSpectra: state.hiddenSpectra ?? [],
          excludeSubtypes: state.excludeSubtypes ?? defaults.excludeSubtypes,
          exNorm: state.exNorm ?? defaults.exNorm,
          chartOptions: state.chartOptions
            ? { ...defaults.chartOptions, ...state.chartOptions }
            : defaults.chartOptions,
          customFilters: state.customFilters ?? {},
          customLasers: state.customLasers ?? {},
          overlapCache: state.overlapCache ?? {},
          _urlInitialized: state._urlInitialized ?? false,
        }),
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
        exNorm: state.exNorm,
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

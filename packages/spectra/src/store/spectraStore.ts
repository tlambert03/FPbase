import { create } from "zustand"
import { createJSONStorage, persist } from "zustand/middleware"
import type { ChartOptions, SpectraStore } from "../types"

// Default chart options
const defaultChartOptions: ChartOptions = {
  showY: true,
  showX: true,
  showGrid: false,
  areaFill: false,
  logScale: false,
  scaleEC: false,
  scaleQY: false,
  extremes: null,
  shareTooltip: true,
  palette: "fpbase",
}

export const useSpectraStore = create<SpectraStore>()(
  persist(
    (set, _get) => ({
      // Initial state
      activeSpectra: [],
      activeOverlaps: [],
      excludeSubtypes: [],
      exNorm: null,
      chartOptions: defaultChartOptions,
      customFilters: {},
      customLasers: {},

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

          return { activeSpectra: newSpectra }
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

      // Subtype management
      setExcludeSubtypes: (subtypes) => set({ excludeSubtypes: subtypes }),

      // Excitation normalization
      setExNorm: (norm) => set({ exNorm: norm }),

      // Chart options
      updateChartOptions: (options) =>
        set((state) => ({
          chartOptions: { ...state.chartOptions, ...options },
        })),

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
    }),
    {
      name: "fpbase-spectra-storage",
      storage: createJSONStorage(() => sessionStorage),
      // Only persist certain fields
      partialize: (state) => ({
        activeSpectra: state.activeSpectra,
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

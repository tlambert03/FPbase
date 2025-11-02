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

// Helper to generate unique selector IDs
let selectorIdCounter = 0
const generateSelectorId = () => `selector-${++selectorIdCounter}`

export const useSpectraStore = create<SpectraStore>()(
  persist(
    (set, _get) => ({
      // Initial state
      activeSpectra: [],
      activeOverlaps: [],
      selectors: [],
      excludeSubtypes: [],
      exNorm: [],
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

      // Selector management
      addSelectors: (selectors) =>
        set((state) => {
          const newSelectors = selectors.map((s) => ({
            ...s,
            id: s.id || generateSelectorId(),
          }))
          return { selectors: [...state.selectors, ...newSelectors] }
        }),

      updateSelector: (selector) =>
        set((state) => ({
          selectors: state.selectors.map((s) => (s.id === selector.id ? { ...s, ...selector } : s)),
        })),

      removeSelector: (id) =>
        set((state) => ({
          selectors: state.selectors.filter((s) => s.id !== id),
        })),

      normalizeCurrent: () =>
        set((state) => {
          const { activeSpectra, selectors } = state

          // Get all unique owners from active spectra
          // This is a simplified version - in reality we'd need to fetch spectrum details
          // to know their owners. For now, we'll just ensure we have enough selectors.

          // Get all owners that have active spectra
          const _ownersWithSpectra = new Set<string>()

          // In the real implementation, we'd iterate through actual spectrum data
          // For now, we'll keep existing selectors that have owners

          // Ensure we have at least one selector
          if (selectors.length === 0 && activeSpectra.length > 0) {
            return {
              selectors: [
                {
                  id: generateSelectorId(),
                  owner: null,
                  category: null,
                },
              ],
            }
          }

          // Remove selectors that don't have owners and aren't the last one
          const hasEmptySelector = selectors.some((s) => !s.owner)
          if (!hasEmptySelector && selectors.length > 0) {
            return {
              selectors: [
                ...selectors,
                {
                  id: generateSelectorId(),
                  owner: null,
                  category: null,
                },
              ],
            }
          }

          return state
        }),

      // Subtype management
      setExcludeSubtypes: (subtypes) => set({ excludeSubtypes: subtypes }),

      // Excitation normalization
      setExNorm: (ids) => set({ exNorm: ids }),

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
          // biome-ignore lint/correctness/noUnusedVariables: removed is needed for destructuring
          const { [id]: removed, ...rest } = state.customFilters
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
          // biome-ignore lint/correctness/noUnusedVariables: removed is needed for destructuring
          const { [id]: removed, ...rest } = state.customLasers
          return { customLasers: rest }
        }),

      // Clear form
      clearForm: (leave = [], appendSpectra = []) =>
        set((state) => {
          let newSpectra = [...state.activeSpectra]
          let newSelectors = [...state.selectors]

          if (leave.length === 0) {
            // Clear everything
            newSpectra = []
            newSelectors = []
          } else {
            // Keep only spectra from specified categories
            // This is simplified - in reality we'd need to check each spectrum's category
            // For now, we'll just keep selectors with matching categories
            newSelectors = newSelectors.filter((s) => s.category && leave.includes(s.category))
          }

          // Add any spectra to append
          if (appendSpectra.length > 0) {
            const toAdd = appendSpectra.filter((id) => !newSpectra.includes(id))
            newSpectra = [...newSpectra, ...toAdd]
          }

          // Ensure at least one selector
          if (newSelectors.length === 0) {
            newSelectors = [
              {
                id: generateSelectorId(),
                owner: null,
                category: null,
              },
            ]
          }

          return {
            activeSpectra: newSpectra,
            selectors: newSelectors,
          }
        }),
    }),
    {
      name: "fpbase-spectra-storage",
      storage: createJSONStorage(() => sessionStorage),
      // Only persist certain fields
      partialize: (state) => ({
        activeSpectra: state.activeSpectra,
        selectors: state.selectors,
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
export const useSelectors = () => useSpectraStore((state) => state.selectors)
export const useExcludeSubtypes = () => useSpectraStore((state) => state.excludeSubtypes)

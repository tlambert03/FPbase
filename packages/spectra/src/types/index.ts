// Core Spectrum Types
export type SpectrumCategory = "P" | "D" | "F" | "L" | "C" | "BP" | "LP" | "SP"
export type SpectrumSubtype = "ex" | "em" | "ab" | "2p" | "bx" | "pd" | "qd"

export interface SpectrumOwner {
  id: string
  slug: string
  name: string
  url?: string
  // Fluorophore-specific fields (for P and D types)
  qy?: number
  extCoeff?: number
  twopPeakgm?: number
  exMax?: number
  emMax?: number
}

export interface Spectrum {
  id: string
  data: [number, number][] // [wavelength, intensity] pairs
  category: SpectrumCategory
  subtype: SpectrumSubtype
  color?: string
  area?: number // Client-side calculated
  owner: SpectrumOwner | null
}

// Custom Spectrum Types (for $cf and $cl)
export interface CustomFilter {
  id: string
  data: [number, number][]
  category: SpectrumCategory
  subtype: SpectrumSubtype
  color?: string
}

export interface CustomLaser {
  id: string
  wavelength: number
  power?: number
  category: SpectrumCategory
}

// Selector Types (for the UI selectors)
export interface Selector {
  id: string
  owner: string | null
  category: SpectrumCategory | null
}

// Chart Options
export interface ChartOptions {
  showY: boolean
  showX: boolean
  showGrid: boolean
  areaFill: boolean
  logScale: boolean
  scaleEC: boolean
  scaleQY: boolean
  extremes: [number, number] | null
  shareTooltip: boolean
  palette: string
}

// Optical Config Types
export interface OpticalConfigFilter {
  id: string
  path: string
  reflects: boolean
  spectrum: {
    id: string
  }
}

export interface OpticalConfigLight {
  id: string
  spectrum: {
    id: string
  }
}

export interface OpticalConfigCamera {
  id: string
  spectrum: {
    id: string
  }
}

export interface OpticalConfig {
  id: number
  name: string
  comments?: string
  microscope: {
    id: string
    name: string
  }
  filters: OpticalConfigFilter[]
  light: OpticalConfigLight | null
  camera: OpticalConfigCamera | null
  laser?: string
}

// Owner Info (for selector options)
export interface OwnerInfo {
  value: string
  label: string
  category: SpectrumCategory
  spectra: Array<{
    id: string
    subtype: SpectrumSubtype
  }>
}

// Store State Types
export interface SpectraState {
  // Active spectra IDs
  activeSpectra: string[]

  // Active overlap IDs
  activeOverlaps: string[]

  // Selectors for UI
  selectors: Selector[]

  // Excluded subtypes
  excludeSubtypes: SpectrumSubtype[]

  // Excitation normalization
  exNorm: string[]

  // Chart options
  chartOptions: ChartOptions

  // Custom spectra (stored in localStorage)
  customFilters: Record<string, CustomFilter>
  customLasers: Record<string, CustomLaser>
}

// Actions
export interface SpectraActions {
  // Spectra management
  setActiveSpectra: (ids: string[]) => void
  updateActiveSpectra: (add?: string[], remove?: string[]) => void

  // Overlaps management
  setActiveOverlaps: (ids: string[]) => void

  // Selector management
  addSelectors: (selectors: Selector[]) => void
  updateSelector: (selector: Selector) => void
  removeSelector: (id: string) => void
  normalizeCurrent: () => void

  // Subtype management
  setExcludeSubtypes: (subtypes: SpectrumSubtype[]) => void

  // Excitation normalization
  setExNorm: (ids: string[]) => void

  // Chart options
  updateChartOptions: (options: Partial<ChartOptions>) => void

  // Custom spectra
  addCustomFilter: (filter: CustomFilter) => void
  removeCustomFilter: (id: string) => void
  addCustomLaser: (laser: CustomLaser) => void
  removeCustomLaser: (id: string) => void

  // Clear form
  clearForm: (leave?: SpectrumCategory[], appendSpectra?: string[]) => void
}

export type SpectraStore = SpectraState & SpectraActions

// GraphQL Response Types
export interface GraphQLResponse<T> {
  data: T
  errors?: Array<{ message: string }>
}

// API Response Types
export interface SpectraListResponse {
  spectra: Array<{
    id: string
    category: SpectrumCategory
    subtype: SpectrumSubtype
    owner: {
      name: string
      slug: string
      url: string
    }
  }>
}

export interface OpticalConfigInfo {
  id: number
  name: string
  comments: string
  microscope: {
    id: string
    name: string
  }
}

export interface OpticalConfigsResponse {
  opticalConfigs: OpticalConfigInfo[]
}

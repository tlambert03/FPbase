// Core Spectrum Types
// Categories returned by GraphQL API (uppercase)
export type SpectrumCategory =
  | "P" // Protein
  | "D" // Dye
  | "F" // Filter
  | "L" // Light
  | "C" // Camera
  | "O" // Overlap (client-side computed)

// Subtypes returned by GraphQL API (uppercase)
// Note: Backend stores these as lowercase, but GraphQL resolver transforms to uppercase
// Graphene-django prefixes "A_" to enum values that aren't valid python identifiers
// (e.g., "2P" starts with a digit), so the API returns "A_2P".
// This is automatically normalized to "2P" at the query boundary (see useSpectraQueries.ts).
// Frontend code should only reference "2P" - "A_2P" is included here only for type safety.
export type SpectrumSubtype =
  | "EX" // Excitation
  | "AB" // Absorption
  | "EM" // Emission
  | "2P" // Two Photon
  | "A_2P" // Two Photon (graphene-django auto-prefixed - normalized to "2P" automatically)
  | "BP" // Bandpass
  | "BX" // Bandpass-Ex
  | "BM" // Bandpass-Em
  | "SP" // Shortpass
  | "LP" // Longpass
  | "BS" // Beamsplitter
  | "QE" // Quantum Efficiency (cameras)
  | "PD" // Power Distribution (lights)
  | "O" // Overlap (client-side computed)

export interface SpectrumOwner {
  id: string
  slug?: string // Optional - not available for all spectrum types (e.g., overlap spectra)
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
  customId?: string // For custom spectra ($cf, $cl)
  data: [number, number][] // [wavelength, intensity] pairs
  category: SpectrumCategory
  subtype: SpectrumSubtype
  color?: string
  area?: number // Client-side calculated
  owner: SpectrumOwner | null
}

// Custom Spectrum Types (for $cf and $cl)

// Parameters for custom filters (stored in state)
export interface CustomFilterParams {
  type: "BP" | "LP" | "SP"
  center: number
  width: number
  transmission: number
}

// Parameters for custom lasers (stored in state)
export interface CustomLaserParams {
  wavelength: number
}

// Legacy types (deprecated - kept for backward compatibility)
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

// Excitation normalization tuple: [wavelength, spectrumId]
export type ExNorm = readonly [wavelength: number | null, spectrumId: string | null] | null

// Store State Types
export interface SpectraState {
  // Active spectra IDs
  activeSpectra: string[]

  // Active overlap IDs
  activeOverlaps: string[]

  // Hidden spectra IDs (temporarily invisible but not removed)
  hiddenSpectra: string[]

  // Excluded subtypes
  excludeSubtypes: SpectrumSubtype[]

  // Excitation normalization: [wavelength, spectrumId]
  exNorm: ExNorm

  // Chart options
  chartOptions: ChartOptions

  // Custom spectra parameters (stored in sessionStorage)
  // Keys are stable IDs like "$cf0", "$cl1"
  customFilters: Record<string, CustomFilterParams>
  customLasers: Record<string, CustomLaserParams>

  // Overlap cache - computed data, NOT persisted (see partialize)
  overlapCache: Record<string, Spectrum>

  // URL initialization tracking (not persisted)
  _urlInitialized: boolean
}

// Actions
export interface SpectraActions {
  // Spectra management
  setActiveSpectra: (ids: string[]) => void
  updateActiveSpectra: (add?: string[], remove?: string[]) => void

  // Overlaps management
  setActiveOverlaps: (ids: string[]) => void
  updateActiveOverlaps: (add?: string[], remove?: string[]) => void

  // Visibility management
  toggleSpectrumVisibility: (id: string) => void
  setHiddenSpectra: (ids: string[]) => void

  // Subtype management
  setExcludeSubtypes: (subtypes: SpectrumSubtype[]) => void

  // Excitation normalization
  setExNorm: (norm: ExNorm) => void

  // Chart options
  updateChartOptions: (options: Partial<ChartOptions>) => void

  // Custom spectra management
  addCustomFilter: (id: string, params: CustomFilterParams) => void
  updateCustomFilter: (id: string, params: Partial<CustomFilterParams>) => void
  removeCustomFilter: (id: string) => void
  addCustomLaser: (id: string, params: CustomLaserParams) => void
  updateCustomLaser: (id: string, params: Partial<CustomLaserParams>) => void
  removeCustomLaser: (id: string) => void

  // Overlap cache management
  setOverlapCache: (id: string, spectrum: Spectrum) => void
  clearOverlapCache: () => void

  // URL initialization tracking
  setUrlInitialized: (value: boolean) => void
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

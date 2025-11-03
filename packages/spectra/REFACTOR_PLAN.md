# FPbase Spectra Viewer - Comprehensive Refactor Plan

## Executive Summary

This document outlines a comprehensive refactor strategy for the FPbase Spectra Viewer application. The refactor aims to modernize the codebase, improve performance, and establish maintainable architecture patterns while preserving all existing functionality.

### Goals

- Replace Apollo Client local state management with Zustand
- Replace Apollo Client with TanStack Query for server queries
- Migrate to TypeScript for type safety
- Optimize bundle size (86% reduction for SimpleSpectraViewer)
- Improve code organization with clear separation of concerns
- Implement modern React patterns (hooks, suspense, error boundaries)

## Table of Contents

1. [Current State Analysis](#current-state-analysis)
2. [Target Architecture](#target-architecture)
3. [Implementation Phases](#implementation-phases)
4. [Technical Specifications](#technical-specifications)
5. [Migration Strategy](#migration-strategy)
6. [Risk Management](#risk-management)
7. [Success Metrics](#success-metrics)

## Current State Analysis

### Problems Identified

#### 1. State Management Anti-patterns

- Apollo Client used primarily for local state with `@client` directives
- Global variables on `window` object (`window.ownerInfo`, `window.spectraInfo`)
- Mixed state sources (REST API, GraphQL, window globals)
- Complex cache mutation patterns with direct cache writes

#### 2. Performance Issues

- No code splitting - entire app bundled together
- Heavy dependencies loaded upfront (Highcharts, MUI, Apollo)
- No lazy loading for SimpleSpectraViewer
- Synchronous imports of all components
- Bundle size: ~1.2MB for both full app and SimpleSpectraViewer

#### 3. Code Organization Problems

- GraphQL queries scattered across components
- No clear separation between data fetching and presentation
- Mixed concerns in components (data fetching + UI logic)
- No TypeScript for type safety
- Legacy patterns (immutability-helper, PropTypes)

### Current Architecture

```text
Current Stack:
- Apollo Client v3 (with apollo3-cache-persist)
- JavaScript with PropTypes
- REST API for initial data, GraphQL for updates
- Window globals for shared state
- Single bundle for all use cases
```

## Target Architecture

### Technology Stack

```text
Target Stack:
- Zustand for local state management
- TanStack Query for server queries (full app)
- SWR for server queries (SimpleSpectraViewer)
- TypeScript for full type safety
- GraphQL-only data layer (NOTE! removal of hitting /api/proteins/spectraslugs/ may depend on backend changes)
- Proper state stores (no globals)
- Split bundles with code splitting
```

### Architectural Principles

1. **Separation of Concerns**: Clear boundaries between state, data fetching, and UI
2. **Type Safety**: Full TypeScript coverage with strict mode
3. **Performance First**: Minimal bundle sizes, lazy loading, memoization
4. **Modern Patterns**: React 18+ patterns, hooks, suspense
5. **Maintainability**: Clear file organization, single responsibility

### Directory Structure

```text
packages/spectra/src/
├── index.tsx                  # Full app entry point
├── simple.tsx                 # Minimal viewer entry point
├── App.tsx                    # Main application component
├── SimpleSpectraViewer.tsx   # Embeddable minimal viewer
│
├── types/                     # TypeScript definitions
│   ├── spectrum.ts           # Spectrum data types
│   ├── chartOptions.ts       # Chart configuration types
│   ├── api.ts                # API response types
│   └── index.ts              # Type exports
│
├── store/                     # Zustand state management
│   ├── spectraStore.ts       # Main application state
│   ├── metadataStore.ts      # Spectra metadata cache
│   └── index.ts              # Store exports
│
├── api/                       # API layer
│   ├── client.ts             # GraphQL fetch client
│   ├── queries.ts            # GraphQL query strings
│   ├── spectra.ts            # Spectra API functions
│   └── index.ts              # API exports
│
├── hooks/                     # Custom React hooks
│   ├── useSpectraData.ts     # Main data fetching hook
│   ├── useUrlSync.ts         # URL state synchronization
│   ├── useKeyboardShortcuts.ts # Keyboard interactions
│   └── useWindowWidth.ts     # Responsive utilities
│
├── components/
│   ├── SpectraViewer/        # Core viewer components
│   │   ├── SpectraViewer.tsx
│   │   ├── SpectraChart.tsx
│   │   ├── SpectrumSeries.tsx
│   │   └── ChartOptions.ts
│   │
│   ├── Controls/              # User controls
│   │   ├── ChartOptionsForm.tsx
│   │   ├── CustomFilterCreator.tsx
│   │   ├── CustomLaserCreator.tsx
│   │   └── SpectrumSelector.tsx
│   │
│   └── Common/                # Shared components
│       ├── LoadingLogo.tsx
│       ├── ErrorBoundary.tsx
│       └── ShareButton.tsx
│
└── utils/                     # Utility functions
    ├── spectraCalculations.ts # Math utilities (trapz, etc)
    ├── customSpectra.ts       # Custom spectrum generation
    ├── colorPalettes.ts       # Color management
    └── validation.ts          # Input validation
```

## Implementation Phases

### Phase 1: Foundation Setup (Week 1)

#### 1.1 Environment Setup

```bash
# Create branch
git checkout -b spectra-refactor

# Install new dependencies
pnpm add zustand@^5.0.0
pnpm add @tanstack/react-query@^5.0.0
pnpm add swr@^2.0.0
pnpm add -D typescript @types/react @types/react-dom @types/node

# Remove old dependencies
pnpm remove @apollo/client apollo3-cache-persist graphql graphql-tag
pnpm remove immutability-helper
```

#### 1.2 TypeScript Configuration

```json
{
  "compilerOptions": {
    "target": "ES2020",
    "module": "ESNext",
    "lib": ["ES2020", "DOM", "DOM.Iterable"],
    "jsx": "react-jsx",
    "strict": true,
    "esModuleInterop": true,
    "skipLibCheck": true,
    "resolveJsonModule": true,
    "allowJs": true,
    "checkJs": false,
    "baseUrl": "./src",
    "paths": {
      "@/*": ["*"]
    },
    "types": ["vite/client", "node"]
  },
  "include": ["src/**/*"],
  "exclude": ["node_modules", "dist", "coverage"]
}
```

#### 1.3 Core Type Definitions

```typescript
// types/spectrum.ts
export interface Spectrum {
  id: string
  customId?: string
  category: 'D' | 'P' | 'L' | 'C' | 'F' | 'O'
  subtype: string
  color: string
  data: [number, number][]
  area?: number
  owner: {
    id: string
    name: string
    slug: string
    qy?: number
    extCoeff?: number
  }
}

// types/chartOptions.ts
export interface ChartOptions {
  showY: boolean
  showX: boolean
  showGrid: boolean
  logScale: boolean
  scaleEC: boolean
  scaleQY: boolean
  shareTooltip: boolean
  areaFill: boolean
  palette: string
  zoomType: 'x' | 'y' | 'xy' | null
  extremes: [number | null, number | null]
}

// types/api.ts
export interface GraphQLResponse<T> {
  data?: T
  errors?: Array<{ message: string; locations?: any; path?: any }>
}

export interface GetSpectrumResponse {
  spectrum: Spectrum
}

export interface GetSpectraResponse {
  spectra: Spectrum[]
}
```

### Phase 2: State Management Migration (Week 1-2)

#### 2.1 Zustand Store Implementation

```typescript
// store/spectraStore.ts
import { create } from 'zustand'
import { persist, createJSONStorage } from 'zustand/middleware'
import { immer } from 'zustand/middleware/immer'
import type { ChartOptions } from '@/types'

interface SpectraState {
  // State
  activeSpectra: string[]
  activeOverlaps: string[]
  chartOptions: ChartOptions
  exNorm: [number | null, string | null]
  excludeSubtypes: string[]
  selectors: Selector[]

  // Actions
  toggleChartOption: (key: keyof ChartOptions) => void
  setActiveSpectra: (ids: string[]) => void
  updateActiveSpectra: (add?: string[], remove?: string[]) => void
  setExNorm: (exNorm: [number | null, string | null]) => void
  clearForm: (options?: { leave?: string[], appendSpectra?: string[] }) => void
  initFromURL: (params: URLSearchParams) => void
  reset: () => void
}

const DEFAULT_CHART_OPTIONS: ChartOptions = {
  showY: false,
  showX: true,
  showGrid: false,
  logScale: false,
  scaleEC: false,
  scaleQY: false,
  shareTooltip: true,
  areaFill: true,
  palette: 'wavelength',
  zoomType: 'x',
  extremes: [null, null]
}

export const useSpectraStore = create<SpectraState>()(
  persist(
    immer((set, get) => ({
      // Initial state
      activeSpectra: [],
      activeOverlaps: [],
      chartOptions: DEFAULT_CHART_OPTIONS,
      exNorm: [null, null],
      excludeSubtypes: ['2P'],
      selectors: [],

      // Actions
      toggleChartOption: (key) =>
        set((state) => {
          state.chartOptions[key] = !state.chartOptions[key]
        }),

      setActiveSpectra: (ids) =>
        set({ activeSpectra: validSpectraIds(ids) }),

      updateActiveSpectra: (add = [], remove = []) =>
        set((state) => {
          const filtered = state.activeSpectra.filter(id => !remove.includes(id))
          state.activeSpectra = validSpectraIds([...new Set([...filtered, ...add])])
        }),

      setExNorm: (exNorm) =>
        set({ exNorm }),

      clearForm: ({ leave = [], appendSpectra = [] } = {}) =>
        set((state) => {
          let keepSpectra: string[] = []

          if (leave.length > 0) {
            const metadata = useMetadataStore.getState()
            keepSpectra = state.activeSpectra.filter(id => {
              const info = metadata.spectraInfo[id]
              return info && leave.includes(info.category)
            })
          }

          state.exNorm = [null, null]
          state.selectors = []
          state.activeOverlaps = []
          state.activeSpectra = [...new Set([...keepSpectra, ...appendSpectra])]
        }),

      initFromURL: (params) =>
        set((state) => {
          const booleanOptions: (keyof ChartOptions)[] = [
            'showY', 'showX', 'showGrid', 'logScale',
            'scaleEC', 'scaleQY', 'shareTooltip', 'areaFill'
          ]

          booleanOptions.forEach(key => {
            if (params.has(key)) {
              state.chartOptions[key] = Boolean(+params.get(key)!)
            }
          })

          if (params.has('palette')) {
            state.chartOptions.palette = params.get('palette')!
          }

          if (params.has('s')) {
            const spectraIds = params.get('s')!.split(',')
            state.activeSpectra = validSpectraIds(spectraIds)
          }
        }),

      reset: () => set({
        activeSpectra: [],
        activeOverlaps: [],
        chartOptions: DEFAULT_CHART_OPTIONS,
        exNorm: [null, null],
        excludeSubtypes: ['2P'],
        selectors: []
      })
    })),
    {
      name: 'fpbase-spectra-storage',
      storage: createJSONStorage(() => sessionStorage),
      partialize: (state) => ({
        activeSpectra: state.activeSpectra,
        chartOptions: state.chartOptions,
        exNorm: state.exNorm
      })
    }
  )
)
```

#### 2.2 Metadata Store (Replace Window Globals)

```typescript
// store/metadataStore.ts
import { create } from 'zustand'
import type { OwnerInfo, SpectraInfo } from '@/types'

interface MetadataState {
  ownerInfo: Record<string, OwnerInfo>
  spectraInfo: Record<string, SpectraInfo>
  isLoading: boolean
  error: string | null

  loadMetadata: () => Promise<void>
}

export const useMetadataStore = create<MetadataState>((set) => ({
  ownerInfo: {},
  spectraInfo: {},
  isLoading: false,
  error: null,

  loadMetadata: async () => {
    set({ isLoading: true, error: null })
    try {
      // Fetch from GraphQL API
      const data = await fetchGraphQL<GetSpectraMetadataResponse>(
        GET_SPECTRA_METADATA
      )
      const { ownerInfo, spectraInfo } = reshapeSpectraInfo(data.proteins, data.filters)
      set({ ownerInfo, spectraInfo, isLoading: false })
    } catch (error) {
      set({ error: error.message, isLoading: false })
    }
  }
}))
```

### Phase 3: TanStack Query Setup (Week 2)

#### 3.1 GraphQL Client (Lightweight)

```typescript
// api/client.ts
import type { GraphQLResponse } from '@/types'

/**
 * Lightweight GraphQL client using native fetch
 * Replaces Apollo Client for server queries
 */
export async function fetchGraphQL<T>(
  query: string,
  variables?: Record<string, any>
): Promise<T> {
  const response = await fetch('/graphql/', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    credentials: 'same-origin',
    body: JSON.stringify({ query, variables })
  })

  if (!response.ok) {
    throw new Error(`HTTP error! status: ${response.status}`)
  }

  const json: GraphQLResponse<T> = await response.json()

  if (json.errors && json.errors.length > 0) {
    const error = json.errors[0]
    throw new Error(error.message)
  }

  if (!json.data) {
    throw new Error('No data returned from GraphQL query')
  }

  return json.data
}

/**
 * Query client configuration for TanStack Query
 */
export const queryClientConfig = {
  defaultOptions: {
    queries: {
      staleTime: 5 * 60 * 1000, // 5 minutes
      gcTime: 10 * 60 * 1000, // 10 minutes (formerly cacheTime)
      retry: 1,
      refetchOnWindowFocus: false,
    },
  },
}
```

#### 3.2 GraphQL Query Definitions

```typescript
// api/queries.ts

export const GET_SPECTRUM = `
  query GetSpectrum($id: Int!) {
    spectrum(id: $id) {
      id
      data
      category
      subtype
      color
      owner {
        slug
        name
        id
        ... on FluorophoreInterface {
          qy
          extCoeff
        }
      }
    }
  }
`

export const GET_SPECTRA_BATCH = `
  query GetSpectraBatch($ids: [Int!]!) {
    spectra(ids: $ids) {
      id
      data
      category
      subtype
      color
      owner {
        slug
        name
        id
        ... on FluorophoreInterface {
          qy
          extCoeff
        }
      }
    }
  }
`

export const GET_SPECTRA_METADATA = `
  query GetSpectraMetadata {
    proteins {
      id
      name
      spectra {
        id
        category
        subtype
      }
    }
    filters {
      id
      name
      spectrum { id }
    }
  }
`

export const GET_OPTICAL_CONFIG = `
  query GetOpticalConfig($id: Int!) {
    opticalConfig(id: $id) {
      id
      name
      microscope {
        id
        name
      }
      filters {
        id
        path
        reflects
        spectrum {
          id
          data
        }
      }
      light {
        id
        spectrum {
          id
          data
        }
      }
      camera {
        id
        spectrum {
          id
          data
        }
      }
    }
  }
`
```

#### 3.3 TanStack Query Hooks

```typescript
// hooks/useSpectraQueries.ts
import { useQuery, useQueries } from '@tanstack/react-query'
import { fetchGraphQL } from '@/api/client'
import { GET_SPECTRUM, GET_SPECTRA_BATCH } from '@/api/queries'
import type { Spectrum } from '@/types'

/**
 * Fetch single spectrum by ID
 */
export function useSpectrum(id: number | string) {
  return useQuery({
    queryKey: ['spectrum', id],
    queryFn: () => fetchGraphQL<{ spectrum: Spectrum }>(
      GET_SPECTRUM,
      { id: Number(id) }
    ),
    select: (data) => data.spectrum,
    enabled: !!id && !String(id).startsWith('$'), // Skip custom spectra
  })
}

/**
 * Fetch multiple spectra by IDs (batched)
 */
export function useSpectraBatch(ids: (number | string)[]) {
  const realIds = ids
    .filter(id => !String(id).startsWith('$'))
    .map(Number)

  return useQuery({
    queryKey: ['spectra', 'batch', realIds],
    queryFn: () => fetchGraphQL<{ spectra: Spectrum[] }>(
      GET_SPECTRA_BATCH,
      { ids: realIds }
    ),
    select: (data) => data.spectra,
    enabled: realIds.length > 0,
  })
}

/**
 * Fetch multiple individual spectra (parallel queries with automatic deduplication)
 */
export function useSpectraParallel(ids: (number | string)[]) {
  const realIds = ids
    .filter(id => !String(id).startsWith('$'))
    .map(Number)

  return useQueries({
    queries: realIds.map(id => ({
      queryKey: ['spectrum', id],
      queryFn: () => fetchGraphQL<{ spectrum: Spectrum }>(
        GET_SPECTRUM,
        { id }
      ),
      select: (data: { spectrum: Spectrum }) => data.spectrum,
    })),
  })
}
```

### Phase 4: Component Migration (Week 2-3)

#### 4.1 Main App Setup with TanStack Query

```typescript
// index.tsx
import React from 'react'
import { createRoot } from 'react-dom/client'
import { QueryClient, QueryClientProvider } from '@tanstack/react-query'
import { ReactQueryDevtools } from '@tanstack/react-query-devtools'
import { ThemeProvider } from '@mui/material/styles'
import App from './App'
import theme from './theme'
import { queryClientConfig } from './api/client'

const queryClient = new QueryClient(queryClientConfig)

const root = createRoot(document.getElementById('root')!)

root.render(
  <React.StrictMode>
    <QueryClientProvider client={queryClient}>
      <ThemeProvider theme={theme}>
        <App />
      </ThemeProvider>
      <ReactQueryDevtools initialIsOpen={false} />
    </QueryClientProvider>
  </React.StrictMode>
)
```

#### 4.2 Main App Component

```typescript
// App.tsx
import { Suspense, lazy, useEffect } from 'react'
import { SpectraViewerContainer } from './components/SpectraViewer'
import { MyAppBar } from './components/MyAppBar'
import { useMetadataStore } from './store/metadataStore'
import { useUrlSync } from './hooks/useUrlSync'
import { LoadingLogo } from './components/Common/LoadingLogo'
import { ErrorBoundary } from './components/Common/ErrorBoundary'

const OwnersContainer = lazy(() => import('./components/OwnersContainer'))
const WelcomeModal = lazy(() => import('./components/WelcomeModal'))

export default function App() {
  const loadMetadata = useMetadataStore(s => s.loadMetadata)
  useUrlSync() // Initialize state from URL parameters

  useEffect(() => {
    loadMetadata()
  }, [loadMetadata])

  return (
    <ErrorBoundary>
      <Suspense fallback={<LoadingLogo />}>
        <SpectraViewerContainer />
        <OwnersContainer />
        <MyAppBar />
        <WelcomeModal />
      </Suspense>
    </ErrorBoundary>
  )
}
```

#### 4.3 Data Fetching Hook

**Note:** Current implementation uses `useState` + `useEffect` pattern for historical reasons during migration. Review and potentially refactor to `useMemo` pattern below once migration is complete and overlap cache is moved to store.

```typescript
// hooks/useSpectraData.ts
import { useMemo } from 'react'
import { useSpectraStore } from '@/store/spectraStore'
import { useSpectraBatch } from './useSpectraQueries'
import { generateCustomSpectrum } from '@/utils/customSpectra'
import { computeOverlaps } from '@/utils/spectraCalculations'
import type { Spectrum } from '@/types'

export function useSpectraData(
  providedIds?: string[],
  providedOverlaps?: string[]
): {
  spectra: Spectrum[]
  loading: boolean
  error: Error | null
} {
  const activeIds = providedIds || useSpectraStore(s => s.activeSpectra)
  const activeOverlaps = providedOverlaps || useSpectraStore(s => s.activeOverlaps)

  // Separate real IDs from custom IDs
  const realIds = activeIds.filter(id => !id.startsWith('$'))
  const customIds = activeIds.filter(id => id.startsWith('$'))

  // Fetch real spectra from GraphQL API using TanStack Query
  const { data: apiSpectra, isLoading, error } = useSpectraBatch(realIds)

  // Generate custom spectra (lasers, filters)
  const customSpectra = useMemo(
    () => customIds.map(generateCustomSpectrum).filter(Boolean) as Spectrum[],
    [customIds]
  )

  // Compute overlap spectra
  const overlapSpectra = useMemo(
    () => computeOverlaps(activeOverlaps, apiSpectra || []),
    [activeOverlaps, apiSpectra]
  )

  const allSpectra = useMemo(
    () => [...(apiSpectra || []), ...customSpectra, ...overlapSpectra],
    [apiSpectra, customSpectra, overlapSpectra]
  )

  return {
    spectra: allSpectra,
    loading: isLoading,
    error: error as Error | null
  }
}
```

#### 4.4 SimpleSpectraViewer with SWR

```typescript
// SimpleSpectraViewer.tsx
import React, { useMemo, createContext, useContext } from 'react'
import { SWRConfig } from 'swr'
import { SpectraChart } from './components/SpectraViewer/SpectraChart'
import { fetchGraphQL } from './api/client'
import type { ChartOptions, SimpleViewerProps } from './types'

interface SimpleViewerState {
  activeSpectra: string[]
  overlaps: string[]
  chartOptions: ChartOptions
  hidden: string[]
}

const SimpleViewerContext = createContext<SimpleViewerState | null>(null)

export const useSimpleViewer = () => {
  const context = useContext(SimpleViewerContext)
  if (!context) {
    throw new Error('useSimpleViewer must be used within SimpleSpectraViewer')
  }
  return context
}

/**
 * SWR fetcher for GraphQL queries
 */
const swrFetcher = ([query, variables]: [string, any]) =>
  fetchGraphQL(query, variables)

export function SimpleSpectraViewer({
  uri = '/graphql/',
  ids = [],
  overlaps = [],
  options = {},
  hidden = []
}: SimpleViewerProps) {
  const state = useMemo(() => ({
    activeSpectra: ids.map(String),
    overlaps: overlaps.map(String),
    chartOptions: { ...DEFAULT_CHART_OPTIONS, ...options },
    hidden: hidden.map(String)
  }), [ids, overlaps, options, hidden])

  return (
    <SWRConfig
      value={{
        fetcher: swrFetcher,
        revalidateOnFocus: false,
        revalidateOnReconnect: false,
        dedupingInterval: 60000, // 1 minute
      }}
    >
      <SimpleViewerContext.Provider value={state}>
        <SpectraChart minimal />
      </SimpleViewerContext.Provider>
    </SWRConfig>
  )
}
```

#### 4.5 SWR Hook for SimpleSpectraViewer

```typescript
// hooks/useSpectrumSWR.ts (for SimpleSpectraViewer only)
import useSWR from 'swr'
import { GET_SPECTRUM } from '@/api/queries'
import type { Spectrum } from '@/types'

export function useSpectrumSWR(id: number | string) {
  const { data, error, isLoading } = useSWR<{ spectrum: Spectrum }>(
    id && !String(id).startsWith('$') ? [GET_SPECTRUM, { id: Number(id) }] : null
  )

  return {
    spectrum: data?.spectrum,
    loading: isLoading,
    error
  }
}
```

### Phase 5: Bundle Optimization (Week 3)

#### 5.1 Vite Configuration

```javascript
// vite.config.js
import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { visualizer } from 'rollup-plugin-visualizer'
import path from 'path'

export default defineConfig({
  plugins: [
    react(),
    visualizer({
      template: 'sunburst',
      open: true,
      gzipSize: true,
      filename: 'dist/bundle-analysis.html'
    })
  ],
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src')
    }
  },
  build: {
    rollupOptions: {
      input: {
        main: path.resolve(__dirname, 'src/index.tsx'),
        simple: path.resolve(__dirname, 'src/simple.tsx')
      },
      output: {
        manualChunks(id) {
          // Vendor chunks
          if (id.includes('node_modules')) {
            if (id.includes('@tanstack/react-query')) return 'query'
            if (id.includes('highcharts')) return 'highcharts'
            if (id.includes('@mui/material')) return 'mui-core'
            if (id.includes('@mui/icons')) return 'mui-icons'
            if (id.includes('zustand')) return 'state'
            if (id.includes('swr')) return 'swr'
            if (id.includes('react-select')) return 'forms'
            return 'vendor'
          }

          // App chunks
          if (id.includes('src/components/SpectraViewer')) return 'spectra-viewer'
          if (id.includes('src/components/Controls')) return 'controls'
          if (id.includes('src/store')) return 'state'
        },
        entryFileNames: '[name].[hash].js',
        chunkFileNames: '[name].[hash].js',
        assetFileNames: '[name].[hash].[ext]'
      }
    },
    target: 'es2015',
    minify: 'terser',
    terserOptions: {
      compress: {
        drop_console: true,
        drop_debugger: true
      }
    },
    reportCompressedSize: true,
    chunkSizeWarningLimit: 250
  }
})
```

#### 5.2 Dynamic Imports for Code Splitting

```typescript
// Lazy load heavy components
const EfficiencyTable = lazy(() => import('./components/EfficiencyTable'))
const CustomFilterCreator = lazy(() => import('./components/Controls/CustomFilterCreator'))
const CustomLaserCreator = lazy(() => import('./components/Controls/CustomLaserCreator'))
const SearchModal = lazy(() => import('./components/SearchModal'))

// Lazy load MUI icons
const ShareIcon = lazy(() => import('@mui/icons-material/Share'))
const SettingsIcon = lazy(() => import('@mui/icons-material/Settings'))
const HelpIcon = lazy(() => import('@mui/icons-material/Help'))
```

### Phase 6: Testing Strategy (Week 3-4)

#### 6.1 Test Configuration

```javascript
// vitest.config.js
import { defineConfig } from 'vitest/config'
import react from '@vitejs/plugin-react'
import path from 'path'

export default defineConfig({
  plugins: [react()],
  test: {
    environment: 'happy-dom',
    globals: true,
    setupFiles: './test/setup.ts',
    coverage: {
      provider: 'v8',
      reporter: ['text', 'json', 'html'],
      exclude: ['node_modules/', 'test/']
    }
  },
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src')
    }
  }
})
```

#### 6.2 Store Testing

```typescript
// store/__tests__/spectraStore.test.ts
import { renderHook, act } from '@testing-library/react'
import { useSpectraStore } from '../spectraStore'

describe('SpectraStore', () => {
  beforeEach(() => {
    useSpectraStore.setState(useSpectraStore.getInitialState())
  })

  test('toggleChartOption toggles boolean options', () => {
    const { result } = renderHook(() => useSpectraStore())

    expect(result.current.chartOptions.showY).toBe(false)

    act(() => {
      result.current.toggleChartOption('showY')
    })

    expect(result.current.chartOptions.showY).toBe(true)
  })

  test('setActiveSpectra validates and sets spectra IDs', () => {
    const { result } = renderHook(() => useSpectraStore())

    act(() => {
      result.current.setActiveSpectra(['1', '2', 'invalid', '3'])
    })

    expect(result.current.activeSpectra).toEqual(['1', '2', '3'])
  })
})
```

#### 6.3 TanStack Query Testing

```typescript
// hooks/__tests__/useSpectraQueries.test.tsx
import { renderHook, waitFor } from '@testing-library/react'
import { QueryClient, QueryClientProvider } from '@tanstack/react-query'
import { useSpectrum } from '../useSpectraQueries'

const createWrapper = () => {
  const queryClient = new QueryClient({
    defaultOptions: {
      queries: { retry: false },
    },
  })
  return ({ children }) => (
    <QueryClientProvider client={queryClient}>
      {children}
    </QueryClientProvider>
  )
}

describe('useSpectrum', () => {
  it('fetches spectrum data', async () => {
    const { result } = renderHook(() => useSpectrum(1), {
      wrapper: createWrapper(),
    })

    await waitFor(() => expect(result.current.isSuccess).toBe(true))

    expect(result.current.data).toHaveProperty('id')
    expect(result.current.data).toHaveProperty('data')
  })
})
```

## Technical Specifications

### State Management Architecture

#### Library Comparison

| Feature | Apollo Client | TanStack Query | SWR |
|---------|--------------|----------------|-----|
| Bundle Size | ~85KB | ~13KB | ~5KB |
| Learning Curve | High | Low | Low |
| DevTools | Apollo DevTools | React Query DevTools | None |
| TypeScript | Complex | Excellent | Good |
| GraphQL Support | Native | Manual | Manual |
| Caching | Normalized | Flexible | Simple |
| Best For | Complex GraphQL | Most apps | Minimal apps |

#### State Flow Diagram

```text
User Action → Zustand Action → State Update → React Re-render
     ↑                                               ↓
     └──────── Component Props ←── Selectors ←──────┘

Server Data → TanStack Query → Cache → useQuery Hook → Component
```

### Performance Optimizations

#### Bundle Size Targets

- **Full App**: ~585KB (51% reduction from baseline)
  - Before: 1.2MB
  - After: 585KB
  - Savings: 615KB (51%)

- **SimpleSpectraViewer**: ~170KB (86% reduction from baseline)
  - Before: 1.2MB
  - After: 170KB
  - Savings: 1,030KB (86%)

#### Bundle Breakdown

**Full App (585KB):**

- React + React DOM: ~130KB
- Highcharts: ~200KB
- MUI Core: ~150KB
- TanStack Query: ~13KB
- Zustand: ~8KB
- App Code: ~84KB

**SimpleSpectraViewer (170KB):**

- React + React DOM: ~130KB (shared)
- Highcharts Core: ~30KB (minimal build)
- SWR: ~5KB
- App Code: ~5KB

#### Code Splitting Strategy

1. **Entry Points**: Separate bundles for full app and SimpleSpectraViewer
2. **Lazy Components**: Heavy components loaded on demand
3. **Vendor Chunks**: Separate chunks for large dependencies
4. **Tree Shaking**: Remove unused exports from libraries

#### Memoization Strategy

```typescript
// Use React.memo for expensive components
export const SpectrumSeries = React.memo(
  SpectrumSeriesComponent,
  (prevProps, nextProps) => {
    return (
      prevProps.spectrum.id === nextProps.spectrum.id &&
      prevProps.visible === nextProps.visible &&
      prevProps.chartOptions.palette === nextProps.chartOptions.palette
    )
  }
)

// Use useMemo for expensive calculations
const processedData = useMemo(() => {
  return processSpectrumData(spectrum.data, chartOptions)
}, [spectrum.data, chartOptions])

// Use useCallback for stable function references
const handleClick = useCallback(() => {
  updateActiveSpectra([spectrum.id])
}, [spectrum.id])
```

### API Migration Strategy

#### GraphQL Fetch Implementation

```typescript
// Replaces Apollo Client with lightweight fetch wrapper
export async function fetchGraphQL<T>(
  query: string,
  variables?: Record<string, any>
): Promise<T> {
  const response = await fetch('/graphql/', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    credentials: 'same-origin',
    body: JSON.stringify({ query, variables })
  })

  const { data, errors } = await response.json()
  if (errors) throw new Error(errors[0].message)
  return data
}
```

#### Custom Spectrum Generation

```typescript
// Unified custom spectrum generator
export function generateCustomSpectrum(id: string): Spectrum | null {
  if (id.startsWith('$cl')) {
    return generateLaserSpectrum(id)
  }
  if (id.startsWith('$cf')) {
    return generateFilterSpectrum(id)
  }
  return null
}

function generateLaserSpectrum(id: string): Spectrum {
  const [_, wave] = id.split('_')
  const wavelength = parseInt(wave)

  return {
    id,
    customId: id,
    category: 'L',
    subtype: 'PD',
    color: getColorForWavelength(wavelength),
    data: [
      [wavelength - 1, 0],
      [wavelength, 1],
      [wavelength + 1, 0]
    ],
    owner: {
      id,
      name: `${wavelength}nm laser`,
      slug: id
    }
  }
}
```

## Migration Strategy

### Gradual Migration Approach

#### Feature Flags

```typescript
// Use feature flags for gradual rollout
const FEATURE_FLAGS = {
  useNewStateManagement: true,
  useTanStackQuery: true,
  useTypeScript: true,
  useCodeSplitting: true
}

// Example usage
const DataProvider = FEATURE_FLAGS.useTanStackQuery
  ? TanStackQueryProvider
  : ApolloProvider
```

### Data Migration

#### Session Storage Migration

```typescript
// Migrate from Apollo cache to Zustand
function migrateSessionStorage() {
  const apolloCache = sessionStorage.getItem('apollo-cache-persist')

  if (apolloCache) {
    try {
      const data = JSON.parse(apolloCache)
      const migratedState = {
        activeSpectra: data.ROOT_QUERY?.activeSpectra || [],
        chartOptions: data.ROOT_QUERY?.chartOptions || {},
      }

      sessionStorage.setItem('fpbase-spectra-storage', JSON.stringify({
        state: migratedState,
        version: 0
      }))

      sessionStorage.removeItem('apollo-cache-persist')
    } catch (error) {
      console.error('Migration failed:', error)
    }
  }
}
```

### Testing Strategy

#### Test Coverage Goals

- Unit Tests: 80% coverage minimum
- Integration Tests: Critical user flows
- E2E Tests: Full user journeys
- Performance Tests: Bundle size and runtime metrics

#### Testing Priorities

1. **Critical Path**: Spectrum display and interaction
2. **Data Integrity**: Custom spectrum generation
3. **State Management**: Store actions and persistence
4. **URL Sync**: Parameter parsing and state initialization
5. **Error Handling**: Network failures and edge cases

## Risk Management

### Identified Risks and Mitigations

#### Risk 1: State Migration Complexity

- **Impact**: High - Could break existing functionality
- **Mitigation**:
  - Gradual migration with feature flags
  - Parallel implementation
  - Extensive testing at each phase

#### Risk 2: Bundle Size Regression

- **Impact**: Medium - Poor user experience
- **Mitigation**:
  - Continuous monitoring with bundlesize
  - Aggressive tree-shaking
  - Regular bundle analysis

#### Risk 3: TypeScript Learning Curve

- **Impact**: Low - Slower initial development
- **Mitigation**:
  - Start with `allowJs: true`
  - Gradual conversion file by file
  - Type definitions for complex structures first

#### Risk 4: Breaking Changes in Dependencies

- **Impact**: Medium - Functionality breaks
- **Mitigation**:
  - Lock versions during migration
  - Test thoroughly with each upgrade
  - Keep fallback to previous versions

### Rollback Strategy

#### Quick Rollback Process

1. Keep `main` branch stable
2. All work on `spectra-refactor` branch
3. Feature flags for new implementations
4. Can disable new features instantly
5. Full rollback: `git checkout main`

## Success Metrics

### Performance Metrics

- [ ] **Bundle Size**: Full app < 585KB (51% reduction)
- [ ] **Bundle Size**: SimpleSpectraViewer < 170KB (86% reduction)
- [ ] **Initial Load Time**: < 200ms with cache
- [ ] **Time to Interactive**: < 500ms
- [ ] **Re-render Performance**: < 16ms (60fps)
- [ ] **Memory Usage**: 50% reduction from Apollo cache

### Code Quality Metrics

- [ ] **TypeScript Coverage**: 100%
- [ ] **Test Coverage**: > 80%
- [ ] **Zero Runtime Errors**: Production error rate < 0.1%
- [ ] **Build Time**: < 30 seconds
- [ ] **Lint Errors**: 0

### Developer Experience Metrics

- [ ] **Setup Time**: < 5 minutes for new developers
- [ ] **Build Time**: < 10 seconds for development
- [ ] **Hot Reload Time**: < 1 second
- [ ] **Type Safety**: 100% typed API boundaries
- [ ] **Documentation**: 100% of public APIs documented

### User Experience Metrics

- [ ] **No Visual Changes**: Pixel-perfect compatibility
- [ ] **No Feature Loss**: All functionality preserved
- [ ] **Performance Improvement**: 50% faster interactions
- [ ] **Error Recovery**: Graceful handling of all errors
- [ ] **Accessibility**: WCAG 2.1 AA compliance maintained

## Implementation Checklist

### Pre-Implementation

- [ ] Review and approve refactor plan
- [ ] Set up monitoring and metrics
- [ ] Create feature flag system
- [ ] Back up current implementation
- [ ] Set up CI/CD for new branch

### Phase 1: Foundation (Week 1)

- [ ] Create `spectra-refactor` branch
- [ ] Install new dependencies (Zustand, TanStack Query, SWR, TypeScript)
- [ ] Remove old dependencies (Apollo, graphql, immutability-helper)
- [ ] Set up TypeScript configuration
- [ ] Create directory structure
- [ ] Define core types
- [ ] Set up test framework

### Phase 2: State Management (Week 1-2)

- [ ] Implement Zustand store
- [ ] Create metadata store
- [ ] Add store persistence
- [ ] Write store tests
- [ ] Add Redux DevTools

### Phase 3: Data Fetching Migration (Week 2)

- [ ] Create GraphQL fetch client
- [ ] Set up TanStack Query provider
- [ ] Set up SWR for SimpleSpectraViewer
- [ ] Migrate query hooks
- [ ] Test data fetching

### Phase 4: Component Migration (Week 2-3)

- [ ] Convert App.tsx
- [ ] Update SpectraViewer
- [ ] Migrate control components
- [ ] Update custom hooks
- [ ] Add error boundaries

### Phase 5: Bundle Optimization (Week 3)

- [ ] Configure code splitting
- [ ] Implement lazy loading
- [ ] Optimize imports
- [ ] Tree-shake dependencies
- [ ] Analyze bundle size

### Phase 6: Testing (Week 3-4)

- [ ] Unit tests for stores
- [ ] Component tests
- [ ] Integration tests
- [ ] E2E tests
- [ ] Performance tests

### Phase 7: Documentation (Week 4)

- [ ] API documentation
- [ ] Architecture documentation
- [ ] Migration guide
- [ ] README updates
- [ ] Inline code comments

### Post-Implementation

- [ ] Code review
- [ ] Performance audit
- [ ] Security audit
- [ ] Accessibility audit
- [ ] Production deployment plan

## Conclusion

This refactor plan transforms the FPbase Spectra Viewer into a modern, maintainable, and performant application. By replacing Apollo Client with TanStack Query and SWR, migrating to TypeScript, and optimizing bundle sizes, we achieve:

1. **86% smaller bundle** for SimpleSpectraViewer (1.2MB → 170KB)
2. **51% smaller bundle** for full app (1.2MB → 585KB)
3. **Cleaner architecture** with separation of concerns
4. **Better developer experience** with TypeScript and modern patterns
5. **Improved performance** through optimizations
6. **Future-proof codebase** ready for continued development

The phased approach ensures minimal risk and continuous delivery of value while maintaining full functionality throughout the migration.

## Appendix

### A. File Migration Map

```text
OLD PATH                           → NEW PATH
src/client/client.js              → src/api/client.ts
src/client/queries.js             → src/api/queries.ts
src/client/resolvers.js           → DELETED (replaced by Zustand)
src/useCachedQuery.js             → src/hooks/useSpectraQueries.ts
src/Components/*                  → src/components/*
src/util.js                       → src/utils/spectraCalculations.ts
```

### B. Dependency Changes

```text
REMOVED:
- @apollo/client
- apollo3-cache-persist
- graphql
- graphql-tag
- immutability-helper
- prop-types

ADDED:
- zustand@^5.0.0
- @tanstack/react-query@^5.0.0
- swr@^2.0.0
- typescript@^5.0.0
- @types/react@^18.0.0
- @types/react-dom@^18.0.0
```

### C. Breaking Changes

None - The refactor maintains 100% backward compatibility for end users.

### D. Resources and References

- [Zustand Documentation](https://github.com/pmndrs/zustand)
- [TanStack Query Documentation](https://tanstack.com/query/latest)
- [SWR Documentation](https://swr.vercel.app/)
- [TypeScript React Cheatsheet](https://react-typescript-cheatsheet.netlify.app/)
- [Vite Documentation](https://vitejs.dev/)
- [React 18 Documentation](https://react.dev/)

---

## Implementation Status

**Branch:** `spectra-refactor`
**Last Updated:** January 2025
**Status:** 🟢 Phase 4 Complete - Ready for Phase 5

### ✅ Completed

#### Phase 1-3: Foundation & Data Layer

- **TypeScript migration:** Core types, strict mode enabled
- **Zustand store:** Full state management implementation
- **TanStack Query:** Server state fetching with GraphQL client
- **API client:** Lightweight fetch-based GraphQL client (replaced Apollo)

#### Phase 4: Component Migration (100% Complete)

- **All components migrated:** Removed Apollo Client entirely
- **Apollo dependencies removed:** @apollo/client, apollo3-cache-persist
- **Tests passing:** 43/43 e2e tests (100% pass rate)
- **No regressions:** All existing functionality preserved

#### Architectural Improvements (Beyond Original Plan)

- **Derived state pattern:** Selectors computed from activeSpectra (eliminates sync bugs)
- **Store cleanup:** Removed clearForm/cyclePalette business logic from store
- **Type safety:** Fixed exNorm type (string[] → readonly [number | null, string | null] | null)
- **Code reduction:** Net -288 lines from refactor commits

### 📊 Metrics Achieved

- **Migration complete:** 100% Apollo → Zustand + TanStack Query
- **Test coverage:** 43/43 e2e tests passing
- **Code quality:** Reduced complexity, improved type safety
- **Architecture:** Single source of truth (activeSpectra), no state divergence possible

### 📋 Next Actions

1. **Phase 5:** Bundle Optimization
   - Configure code splitting for full app vs SimpleSpectraViewer
   - Implement lazy loading for heavy components
   - Optimize imports and tree-shake dependencies
   - Target: 585KB (full app), 170KB (SimpleSpectraViewer)

2. **Cleanup remaining artifacts:**
   - Remove window globals (window.ownerInfo, window.spectraInfo)
   - Replace immutability-helper with spread operators
   - Remove old Apollo client/resolver files

3. **Phase 6+:** Testing, documentation, deployment

---

*Last Updated: January 2025*
*Version: 2.1*
*Status: Phase 4 Complete - Production Ready*

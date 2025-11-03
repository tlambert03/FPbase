import { useMemo } from "react"
import "./index.css"
import { StyledEngineProvider, ThemeProvider } from "@mui/material/styles"
import { QueryClient, QueryClientProvider } from "@tanstack/react-query"
import { ReactQueryDevtools } from "@tanstack/react-query-devtools"
import App from "./App"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import theme from "./Components/theme"
import { queryClientConfig } from "./hooks/useSpectraQueries"
import { useSpectraStore } from "./store/spectraStore"
import { syncURLToStore } from "./store/urlSync"
import type { ChartOptions } from "./types"

const AppWrapper = () => {
  const queryClient = useMemo(() => new QueryClient(queryClientConfig), [])

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <QueryClientProvider client={queryClient}>
          <App />
          <ReactQueryDevtools initialIsOpen={false} buttonPosition="top-right" />
        </QueryClientProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

export default AppWrapper

interface SimpleSpectraViewerProps {
  ids?: string[] | number[]
  overlaps?: string[] | number[]
  options?: Partial<ChartOptions>
  hidden?: string[] | number[]
}

const SimpleSpectraViewer = ({ ids, overlaps, options, hidden }: SimpleSpectraViewerProps) => {
  const queryClient = useMemo(() => new QueryClient(queryClientConfig), [])

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <QueryClientProvider client={queryClient}>
          <Inner ids={ids} overlaps={overlaps} options={options} hidden={hidden} />
          <ReactQueryDevtools initialIsOpen={false} buttonPosition="top-right" />
        </QueryClientProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

interface InnerProps {
  ids?: string[] | number[]
  overlaps?: string[] | number[]
  options?: Partial<ChartOptions>
  hidden?: string[] | number[]
}

const Inner = ({ ids, overlaps, options, hidden }: InnerProps) => {
  // Normalize inputs
  let normalizedIds = Array.isArray(ids) ? ids.map(String) : []
  const normalizedOverlaps = Array.isArray(overlaps) ? overlaps.map(String) : []
  const normalizedHidden = Array.isArray(hidden) ? hidden.map(String) : []
  let normalizedOptions = options || {}

  // If no IDs provided, try to parse from URL
  if (normalizedIds.length === 0) {
    const store = useSpectraStore.getState()
    syncURLToStore(store)
    normalizedIds = store.activeSpectra
    normalizedOptions = { ...store.chartOptions, ...normalizedOptions }
  }

  return (
    <SpectraViewerContainer
      provideSpectra={normalizedIds}
      provideOverlaps={normalizedOverlaps}
      provideOptions={normalizedOptions}
      provideHidden={normalizedHidden}
    />
  )
}

export { SimpleSpectraViewer }

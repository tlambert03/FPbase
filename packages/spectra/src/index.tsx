import { lazy, Suspense, useEffect, useMemo } from "react"
import "./index.css"
import { StyledEngineProvider, ThemeProvider } from "@mui/material/styles"
import { QueryClient, QueryClientProvider } from "@tanstack/react-query"
import App from "./App"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import theme from "./Components/theme"
import { queryClientConfig } from "./hooks/useSpectraQueries"
import { useSpectraStore } from "./store/spectraStore"
import type { SpectraState } from "./types"
import { parseURLParams } from "./utils/urlParams"

// Lazy load devtools only in development - completely excluded from production bundle
const ReactQueryDevtools =
  process.env.NODE_ENV === "development"
    ? lazy(() =>
        import("@tanstack/react-query-devtools").then((mod) => ({
          default: mod.ReactQueryDevtools,
        }))
      )
    : null

const AppWrapper = () => {
  const queryClient = useMemo(() => new QueryClient(queryClientConfig), [])

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <QueryClientProvider client={queryClient}>
          <App />
          {process.env.NODE_ENV === "development" && ReactQueryDevtools && (
            <Suspense fallback={null}>
              <ReactQueryDevtools initialIsOpen={false} buttonPosition="top-right" />
            </Suspense>
          )}
        </QueryClientProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

export default AppWrapper

/**
 * SimpleSpectraViewer - A simplified, view-only version of the spectra viewer.
 *
 * Usage:
 * - `<SimpleSpectraViewer />` or `<SimpleSpectraViewer fromUrl />` - loads from URL (default)
 * - `<SimpleSpectraViewer state={{activeSpectra: ['1', '2'], ...}} />` - uses provided state
 *
 * This component is used for embedding (e.g., spectra_graph.html) and has no UI controls.
 * Must provide EITHER state OR fromUrl, but not both.
 *
 * Implementation: Populates the Zustand store on mount, then SpectraViewerContainer reads from it.
 * This matches the pattern used by App.tsx and simplifies state management.
 */
const SimpleSpectraViewer = ({
  state,
  fromUrl = true,
}: {
  state?: Partial<SpectraState>
  fromUrl?: boolean
}) => {
  const queryClient = useMemo(() => new QueryClient(queryClientConfig), [])

  // Populate the store on mount
  useEffect(() => {
    if (state && fromUrl) {
      throw new Error("SimpleSpectraViewer: Cannot provide both 'state' and 'fromUrl=true'")
    }

    const store = useSpectraStore.getState()

    if (state) {
      // Use provided state
      store.replace(state)
    } else if (fromUrl) {
      // Load from URL parameters
      const urlState = parseURLParams(window.location.search)
      if (Object.keys(urlState).length > 0) {
        store.replace(urlState)
      }
      // If no URL params, keep existing store state (from sessionStorage)
    }
    // If both state and fromUrl are false/undefined, keep existing store state
  }, [state, fromUrl])

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <QueryClientProvider client={queryClient}>
          <SpectraViewerContainer />
          {/* Devtools disabled for SimpleSpectraViewer to prevent unwanted padding-bottom */}
        </QueryClientProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

export { SimpleSpectraViewer }

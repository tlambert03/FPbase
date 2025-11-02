import { useMemo } from "react"
import "./index.css"
import { ApolloProvider } from "@apollo/client"
import { StyledEngineProvider, ThemeProvider } from "@mui/material/styles"
import { QueryClient, QueryClientProvider } from "@tanstack/react-query"
import { ReactQueryDevtools } from "@tanstack/react-query-devtools"
import App from "./App"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import theme from "./Components/theme"
import intializeClient from "./client/client"
import { queryClientConfig } from "./hooks/useSpectraQueries"
import { useSpectraStore } from "./store/spectraStore"
import { syncURLToStore } from "./store/urlSync"
import type { ChartOptions } from "./types"

// biome-ignore lint/correctness/noUnusedFunctionParameters: uri kept for API compatibility
const AppWrapper = ({ uri = "/graphql/" }) => {
  const queryClient = useMemo(() => new QueryClient(queryClientConfig), [])
  const apolloClient = useMemo(() => intializeClient({ uri }), [uri])

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <ApolloProvider client={apolloClient}>
          <QueryClientProvider client={queryClient}>
            <App />
            <ReactQueryDevtools initialIsOpen={false} />
          </QueryClientProvider>
        </ApolloProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

export default AppWrapper

interface SimpleSpectraViewerProps {
  uri?: string
  ids?: string[] | number[]
  overlaps?: string[] | number[]
  options?: Partial<ChartOptions>
  hidden?: string[] | number[]
}

// biome-ignore lint/correctness/noUnusedFunctionParameters: uri kept for API compatibility
const SimpleSpectraViewer = ({
  uri = "/graphql/",
  ids,
  overlaps,
  options,
  hidden,
}: SimpleSpectraViewerProps) => {
  const queryClient = useMemo(() => new QueryClient(queryClientConfig), [])
  const apolloClient = useMemo(() => intializeClient({ uri }), [uri])

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <ApolloProvider client={apolloClient}>
          <QueryClientProvider client={queryClient}>
            <Inner ids={ids} overlaps={overlaps} options={options} hidden={hidden} />
          </QueryClientProvider>
        </ApolloProvider>
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

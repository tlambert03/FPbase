import { useMemo } from "react"
import "./index.css"
import { StyledEngineProvider, ThemeProvider } from "@mui/material/styles"
import { QueryClient, QueryClientProvider } from "@tanstack/react-query"
import { ReactQueryDevtools } from "@tanstack/react-query-devtools"
import App from "./App"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import theme from "./Components/theme"
import { defaults } from "./defaults"
import { queryClientConfig } from "./hooks/useSpectraQueries"
import { parseURL } from "./util"

const AppWrapper = ({ uri = "/graphql/" }) => {
  const queryClient = useMemo(() => new QueryClient(queryClientConfig), [])
  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <QueryClientProvider client={queryClient}>
          <App />
          <ReactQueryDevtools initialIsOpen={false} />
        </QueryClientProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

export default AppWrapper

const SimpleSpectraViewer = ({ uri = "/graphql/", ids, overlaps, options, hidden }) => {
  const queryClient = useMemo(() => new QueryClient(queryClientConfig), [])

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <QueryClientProvider client={queryClient}>
          <Inner ids={ids} overlaps={overlaps} options={options} hidden={hidden} />
          <ReactQueryDevtools initialIsOpen={false} />
        </QueryClientProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

const Inner = ({ ids, overlaps, options, hidden }) => {
  // Normalize null values to empty arrays and options to an empty object if null
  ids = Array.isArray(ids) ? ids : []
  overlaps = Array.isArray(overlaps) ? overlaps : []
  hidden = Array.isArray(hidden) ? hidden : []
  options = options || {}

  if (ids.length === 0) {
    const url_data = parseURL()
    console.log(url_data)
    ids = ids.length > 0 ? ids : url_data.activeSpectra || []
    options = Object.keys(options).length > 0 ? options : url_data.chartOptions || {}
  }

  return (
    <SpectraViewerContainer
      provideSpectra={ids.map(String)}
      provideOverlaps={overlaps}
      provideOptions={{ ...defaults.chartOptions, ...options }}
      provideHidden={hidden.map(String)}
    />
  )
}

export { SimpleSpectraViewer }

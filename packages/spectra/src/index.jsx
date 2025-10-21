import React, { useMemo } from "react"
import "./index.css"
import { ApolloProvider } from "@apollo/client"
import { ThemeProvider, StyledEngineProvider } from "@mui/material/styles"
import App from "./App"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import initializeClient from "./client/client"
import { parseURL } from "./client/client"
import { defaults } from "./client/resolvers"
import theme from "./Components/theme"

/* eslint-disable */
// if (process.env.NODE_ENV === 'development') {
//   const whyDidYouRender = require('@welldone-software/why-did-you-render');
//   whyDidYouRender(React, {
//     trackAllPureComponents: true,
//     logOnDifferentValues: true,
//     collapseGroups: true,
//   })
// }
/* eslint-enable */

const AppWrapper = ({ uri = "/graphql/" }) => {
  const client = useMemo(() => initializeClient({ uri }), [uri])
  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <ApolloProvider client={client}>
          <App />
        </ApolloProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

export default AppWrapper

const SimpleSpectraViewer = ({ uri = "/graphql/", ids, overlaps, options, hidden }) => {
  const client = useMemo(() => initializeClient({ uri }), [uri])

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <ApolloProvider client={client}>
          <Inner ids={ids} overlaps={overlaps} options={options} hidden={hidden} />
        </ApolloProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

const Inner = ({ ids, overlaps, options, hidden }) => {
  // Normalize null values to empty arrays and options to an empty object if null
  ids = Array.isArray(ids) ? ids : [];
  overlaps = Array.isArray(overlaps) ? overlaps : [];
  hidden = Array.isArray(hidden) ? hidden : [];
  options = options || {};

  if (ids.length === 0) {
    const url_data = parseURL();
    console.log(url_data);
    ids = ids.length > 0 ? ids : url_data.activeSpectra || [];
    options = Object.keys(options).length > 0 ? options : url_data.chartOptions || {};
  }

  return (
    <SpectraViewerContainer
      provideSpectra={ids.map(String)}
      provideOverlaps={overlaps}
      provideOptions={{ ...defaults.chartOptions, ...options }}
      provideHidden={hidden.map(String)}
    />
  );
};


export { SimpleSpectraViewer }

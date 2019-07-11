import { ApolloClient } from "apollo-client"
import {
  InMemoryCache,
  IntrospectionFragmentMatcher,
} from "apollo-cache-inmemory"
import { HttpLink } from "apollo-link-http"
import { onError } from "apollo-link-error"
import { ApolloLink } from "apollo-link"
import { persistCache } from "apollo-cache-persist"
import qs from "qs"
import { defaults, resolvers, validSpectraIds } from "./resolvers"
import typeDefs from "./schema"

import introspectionQueryResultData from "../fragmentTypes.json"
import { decoder } from "../util"
import { GET_CHART_OPTIONS } from "./queries"
import "unfetch/polyfill/index"

function intializeClient({ uri, storage }) {
  const fragmentMatcher = new IntrospectionFragmentMatcher({
    introspectionQueryResultData,
  })

  const cache = new InMemoryCache({ fragmentMatcher })

  const link = ApolloLink.from([
    onError(({ graphQLErrors, networkError }) => {
      if (graphQLErrors)
        graphQLErrors.map(({ message, locations, path }) =>
          console.log(
            `[GraphQL error]: Message: ${message}, Location: ${locations}, Path: ${path}`
          )
        )
      if (networkError) console.log(`[Network error]: ${networkError}`)
    }),
    new HttpLink({
      uri: uri || "https://www.fpbase.org/graphql/",
      credentials: "same-origin",
    }),
  ])

  const client = new ApolloClient({
    link,
    cache,
    typeDefs,
    resolvers,
  })

  // Populate from localstorage?
  const setupLocalStorage = async () => {
    cache.writeData({ data: defaults })
    await persistCache({
      cache,
      storage: storage || window.sessionStorage,
      debounce: 400,
    })
    cache.writeData({ data: { activeOverlaps: [] } })
  }

  function parseURL() {
    const url = qs.parse(window.location.search.replace(/^\?/, ""), {
      decoder,
    })
    if (Object.keys(url).length === 0 && url.constructor === Object) return

    const data = cache.readQuery({ query: GET_CHART_OPTIONS })
    data.selectors = []
    const booleanOptions = Object.keys(defaults.chartOptions).filter(
      key => typeof defaults.chartOptions[key] === "boolean"
    )

    const extremes = [null, null]
    const exNorm = [null, null]
    Object.keys(url).forEach(key => {
      if (booleanOptions.includes(key)) {
        data.chartOptions[key] = Boolean(+url[key])
      }
      if (key === "xMin") extremes[0] = +url[key]
      if (key === "xMax") extremes[1] = +url[key]
      if (key === "normWave") exNorm[0] = url[key]
      if (key === "normID") exNorm[1] = url[key]
      if (["s", "activeSpectra"].includes(key)) {
        let active = url[key]
        if (!Array.isArray(active)) active = active.split(",")
        data.activeSpectra = validSpectraIds(active)
      }
    })
    if (extremes.some(i => i)) data.chartOptions.extremes = extremes
    if (exNorm.some(i => i)) data.exNorm = exNorm
    cache.writeData({ data })
  }

  setupLocalStorage().then(parseURL)

  return client
}

export default intializeClient

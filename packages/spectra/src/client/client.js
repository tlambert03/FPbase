import { ApolloClient, InMemoryCache, HttpLink, ApolloLink, from } from "@apollo/client";
import { onError } from "@apollo/client/link/error";
import { persistCache } from "apollo3-cache-persist";
import qs from "qs";
import { defaults, resolvers, validSpectraIds } from "./resolvers";
import typeDefs from "./schema";

import introspectionQueryResultData from "../fragmentTypes.json";
import { decoder } from "../util";
import { GET_CHART_OPTIONS } from "./queries";
import "unfetch/polyfill/index";
import PALETTES from "../palettes";

export function parseURL(data) {
  data = data || defaults;
  const url = qs.parse(window.location.search.replace(/^\?/, ""), {
    decoder,
  });
  if (Object.keys(url).length === 0 && url.constructor === Object) return data;

  data.selectors = [];
  const booleanOptions = Object.keys(defaults.chartOptions).filter(
    (key) => typeof defaults.chartOptions[key] === "boolean"
  );

  const extremes = [null, null];
  const exNorm = [null, null];
  Object.keys(url).forEach((key) => {
    if (booleanOptions.includes(key)) {
      data.chartOptions[key] = Boolean(+url[key]);
    }
    if (key === "palette" && url[key] in PALETTES) {
      data.chartOptions.palette = url[key];
    }
    if (key === "zoomType") {
      if (["x", "y", "xy"].includes(url[key])) {
        data.chartOptions.zoomType = url[key];
      } else {
        data.chartOptions.zoomType = null;
      }
    }
    if (key === "xMin") extremes[0] = +url[key];
    if (key === "xMax") extremes[1] = +url[key];
    if (key === "normWave") exNorm[0] = url[key];
    if (key === "normID") exNorm[1] = url[key];
    if (["s", "activeSpectra"].includes(key)) {
      let active = url[key];
      if (!Array.isArray(active)) active = active.split(",");
      data.activeSpectra = validSpectraIds(active);
    }
  });
  if (extremes.some((i) => i)) data.chartOptions.extremes = extremes;
  if (exNorm.some((i) => i)) data.exNorm = exNorm;
  return data;
}

function intializeClient({ uri, storage }) {
  const cache = new InMemoryCache({
    possibleTypes: introspectionQueryResultData.possibleTypes,
  });

  const link = from([
    onError(({ graphQLErrors, networkError }) => {
      if (graphQLErrors)
        graphQLErrors.map(({ message, locations, path }) =>
          console.log(
            `[GraphQL error]: Message: ${message}, Location: ${locations}, Path: ${path}`
          )
        );
      if (networkError) console.log(`[Network error]: ${networkError}`);
    }),
    new HttpLink({
      uri: uri || "https://www.fpbase.org/graphql/",
      credentials: "same-origin",
    }),
  ]);

  const client = new ApolloClient({
    link,
    cache,
    typeDefs,
    resolvers,
  });

  // Populate from localstorage?
  const setupLocalStorage = async () => {
    try {
      await persistCache({
        cache,
        storage: storage || window.sessionStorage,
        debounce: 400,
        maxSize: 1048576, // 1MB limit to prevent quota exceeded errors
      });
    } catch (error) {
      // If persistence fails (quota exceeded, etc), just continue without it
      console.warn('Cache persistence disabled:', error.message);
    }

    // After restoring from cache, ensure defaults are set if not present
    try {
      const cached = cache.readQuery({ query: GET_CHART_OPTIONS });
      if (!cached || !cached.chartOptions) {
        // No cached data, initialize with defaults
        cache.writeQuery({
          query: GET_CHART_OPTIONS,
          data: {
            chartOptions: defaults.chartOptions,
          },
        });
      }
    } catch (e) {
      // Query failed (no data in cache), initialize with defaults
      cache.writeQuery({
        query: GET_CHART_OPTIONS,
        data: {
          chartOptions: defaults.chartOptions,
        },
      });
    }
  };

  function _parseURL() {
    try {
      const cached = cache.readQuery({ query: GET_CHART_OPTIONS });
      if (!cached) return;

      // parseURL expects the full defaults structure, so we need to merge with defaults
      const fullData = { ...defaults, ...cached };
      const parsedData = parseURL(fullData);
      // Write back only the chartOptions part
      cache.writeQuery({
        query: GET_CHART_OPTIONS,
        data: {
          chartOptions: parsedData.chartOptions,
        },
      });
    } catch (error) {
      console.warn('Failed to parse URL parameters:', error);
    }
  }

  setupLocalStorage().then(_parseURL);

  return client;
}

export default intializeClient;

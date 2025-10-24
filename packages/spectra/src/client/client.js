import { ApolloClient, InMemoryCache, HttpLink, ApolloLink, from } from "@apollo/client";
import { onError } from "@apollo/client/link/error";
import { persistCache } from "apollo3-cache-persist";
import qs from "qs";
import { defaults, resolvers, validSpectraIds } from "./resolvers";
import typeDefs from "./schema";

import introspectionQueryResultData from "../fragmentTypes.json";
import { decoder } from "../util";
import {
  GET_CHART_OPTIONS,
  GET_ACTIVE_SPECTRA,
  GET_ACTIVE_OVERLAPS,
  GET_OWNER_OPTIONS,
  GET_EX_NORM,
  GET_SELECTORS,
} from "./queries";
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
  // Transform Apollo 2.x introspection format to Apollo 3.x possibleTypes format
  // Apollo 2.x format: {__schema: {types: [{name: "Interface", possibleTypes: [{name: "Type"}]}]}}
  // Apollo 3.x format: {Interface: ["Type1", "Type2"]}
  const possibleTypes = introspectionQueryResultData.__schema.types.reduce(
    (acc, type) => {
      if (type.possibleTypes) {
        acc[type.name] = type.possibleTypes.map(t => t.name);
      }
      return acc;
    },
    {}
  );

  const cache = new InMemoryCache({
    possibleTypes,
    // Note: canonizeResults option removed - deprecated in Apollo 3.14+
    // freezeResults: false removed - components should clone data before mutation
    // Freezing cache results prevents bugs from shared state mutations
    typePolicies: {
      Query: {
        fields: {
          chartOptions: {
            merge(existing, incoming) {
              // Merge chartOptions object properly
              return { ...existing, ...incoming };
            },
          },
        },
      },
    },
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
    // typeDefs removed - not supported in Apollo Client v3
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

    // Initialize all cache fields with defaults if not present
    const initializeField = (query, data) => {
      try {
        const cached = cache.readQuery({ query });
        if (!cached) {
          cache.writeQuery({ query, data });
        }
      } catch (e) {
        // Query failed (no data in cache), initialize with defaults
        cache.writeQuery({ query, data });
      }
    };

    initializeField(GET_CHART_OPTIONS, { chartOptions: defaults.chartOptions });
    initializeField(GET_ACTIVE_SPECTRA, { activeSpectra: defaults.activeSpectra });
    initializeField(GET_ACTIVE_OVERLAPS, { activeOverlaps: defaults.activeOverlaps });
    initializeField(GET_OWNER_OPTIONS, { excludeSubtypes: defaults.excludeSubtypes });
    initializeField(GET_EX_NORM, { exNorm: defaults.exNorm });
    initializeField(GET_SELECTORS, { selectors: defaults.selectors });
  };

  function _parseURL() {
    try {
      // Read current chartOptions from cache (restored from sessionStorage)
      // Note: readQuery returns null if the cache is empty, so we need to handle that
      const cachedData = cache.readQuery({ query: GET_CHART_OPTIONS });
      const chartOptions = cachedData?.chartOptions;

      // Create a mutable clone of chartOptions for parseURL to mutate
      // In Apollo v2, objects from cache were mutable; in v3 they're frozen
      const data = {
        chartOptions: chartOptions ? { ...chartOptions } : { ...defaults.chartOptions },
      };

      // parseURL will mutate data and add activeSpectra, exNorm, selectors based on URL params
      // This mimics the original Apollo v2 behavior where parseURL would add fields to the object
      const parsedData = parseURL(data);

      // Write all fields back to cache, mimicking Apollo v2's cache.writeData({ data })
      // which would merge all fields at once
      cache.writeQuery({
        query: GET_CHART_OPTIONS,
        data: { chartOptions: parsedData.chartOptions },
      });

      // Only write activeSpectra if parseURL set it (i.e., if URL had "s" param)
      // Otherwise leave whatever was in sessionStorage
      if (parsedData.activeSpectra !== undefined) {
        cache.writeQuery({
          query: GET_ACTIVE_SPECTRA,
          data: { activeSpectra: parsedData.activeSpectra },
        });
      }

      if (parsedData.exNorm !== undefined) {
        cache.writeQuery({
          query: GET_EX_NORM,
          data: { exNorm: parsedData.exNorm },
        });
      }

      // parseURL always sets selectors to []
      if (parsedData.selectors !== undefined) {
        cache.writeQuery({
          query: GET_SELECTORS,
          data: { selectors: parsedData.selectors },
        });
      }
    } catch (error) {
      console.warn('Failed to parse URL parameters:', error);
    }
  }

  setupLocalStorage().then(_parseURL);

  return client;
}

export default intializeClient;

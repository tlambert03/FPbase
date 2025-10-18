import ApolloClient, { InMemoryCache } from "apollo-boost"
import { persistCache } from "apollo-cache-persist"
import { defaults, resolvers } from "./resolvers"
import { typeDefs } from "./schema"

const cache = new InMemoryCache()
const client = new ApolloClient({
  cache: cache,
  uri: "http://localhost:8000/graphql/",
  clientState: {
    defaults,
    typeDefs,
    resolvers
  }
})
window.client = client

// Populate from localstorage?
const setupLocalStorage = async () => {
  try {
    await persistCache({
      cache,
      storage: window.sessionStorage,
      maxSize: 1048576, // 1MB limit to prevent quota exceeded errors
    })
    // After restoring from cache, reset selectors to empty
    // They'll be regenerated from activeSpectra by NORMALIZE_CURRENT
    cache.writeData({ data: { selectors: [] } });
  } catch (error) {
    // If persistence fails (quota exceeded, etc), just continue without it
    console.warn('Cache persistence disabled:', error.message);
  }
}
setupLocalStorage()

export default client

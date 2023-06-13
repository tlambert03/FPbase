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
  await persistCache({
    cache,
    storage: window.sessionStorage
  })
}
setupLocalStorage()

export default client

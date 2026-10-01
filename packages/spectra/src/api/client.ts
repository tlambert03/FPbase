import type { GraphQLResponse } from "../types"

// longer queries are sent as POST (servers and proxies limit URL length)
const MAX_GET_URL_LENGTH = 2000

/**
 * Lightweight GraphQL client using native fetch
 * Replaces Apollo Client for server queries
 *
 * Queries are sent as GET, which the CDN can cache (all of ours are reads).
 */
export async function fetchGraphQL<T>(
  query: string,
  variables?: Record<string, unknown>
): Promise<T> {
  const params = new URLSearchParams({ query })
  if (variables) params.set("variables", JSON.stringify(variables))
  const url = `/graphql/?${params}`
  const response =
    url.length <= MAX_GET_URL_LENGTH
      ? // (no JSON Content-Type on a GET: the server would try to parse the empty body)
        await fetch(url, { credentials: "same-origin" })
      : await fetch("/graphql/", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          credentials: "same-origin",
          body: JSON.stringify({ query, variables }),
        })

  if (!response.ok) {
    throw new Error(`GraphQL request failed: ${response.status} ${response.statusText}`)
  }

  const json: GraphQLResponse<T> = await response.json()

  if (json.errors?.length) {
    throw new Error(json.errors[0]!.message)
  }

  if (!json.data) {
    throw new Error("GraphQL response missing data")
  }

  return json.data
}

/**
 * Fetch from REST API endpoints
 */
export async function fetchAPI<T>(endpoint: string): Promise<T> {
  const response = await fetch(endpoint, {
    method: "GET",
    credentials: "same-origin",
  })

  if (!response.ok) {
    throw new Error(`API request failed: ${response.status} ${response.statusText}`)
  }

  return response.json()
}

import type { GraphQLResponse } from "../types"

export interface GraphQLOptions {
  variables?: Record<string, unknown>
  method?: "GET" | "POST"
}

/**
 * Lightweight GraphQL client using native fetch
 * Replaces Apollo Client for server queries
 */
export async function fetchGraphQL<T>(query: string, options: GraphQLOptions = {}): Promise<T> {
  const { variables, method = "POST" } = options
  const minified_query = query.replace(/\s+/g, " ").trim()

  let url = "/graphql/"
  const fetchOptions: RequestInit = {
    method,
    credentials: "same-origin",
  }

  if (method === "GET") {
    // Minify and encode query in URL
    url = `/graphql/?query=${encodeURIComponent(minified_query)}`
  } else {
    // POST with JSON body
    fetchOptions.headers = { "Content-Type": "application/json" }
    fetchOptions.body = JSON.stringify({ query: minified_query, variables })
  }

  const response = await fetch(url, fetchOptions)

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

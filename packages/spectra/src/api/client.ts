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

  // Extract operation name from query (e.g., "query SpectraList {" -> "SpectraList")
  const operationMatch = minified_query.match(/^(query|mutation)\s+(\w+)/)
  const operationName = operationMatch?.[2]

  let url = "/graphql/"
  const fetchOptions: RequestInit = {
    method,
    credentials: "same-origin",
  }

  if (method === "GET") {
    // Minify and encode query in URL
    fetchOptions.headers = { "Content-Type": "application/json" }
    const params = new URLSearchParams({ query: minified_query })
    if (operationName) params.set("operationName", operationName)
    url = `/graphql/?${params.toString()}`
  } else {
    // POST with JSON body
    fetchOptions.headers = { "Content-Type": "application/json" }
    fetchOptions.body = JSON.stringify({
      query: minified_query,
      variables,
      ...(operationName && { operationName }),
    })
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

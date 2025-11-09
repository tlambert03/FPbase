import type { GraphQLResponse } from "../types"

/**
 * Normalize GraphQL query by removing unnecessary whitespace.
 *
 * This:
 * - Reduces URL length for GET requests (less encoded whitespace)
 * - Ensures consistent cache keys (same query, different formatting = same URL)
 * - Safe for queries without string literals with significant whitespace
 */
function normalizeGraphQLQuery(query: string): string {
  return (
    query
      .trim()
      // Replace all whitespace sequences (including newlines) with single space
      .replace(/\s+/g, " ")
      // Remove spaces around GraphQL punctuation
      .replace(/\s*([{}(),:])\s*/g, "$1")
  )
}

/**
 * Lightweight GraphQL client using native fetch
 *
 * Uses GET requests for queries (enables browser caching with ETags)
 * and POST for mutations.
 */
export async function fetchGraphQL<T>(
  query: string,
  variables?: Record<string, unknown>,
  options?: { method?: "GET" | "POST" }
): Promise<T> {
  // Default to GET for queries (enables automatic browser caching)
  // Use POST for mutations or when explicitly requested
  const method = options?.method ?? "GET"
  const isMutation = query.trim().startsWith("mutation")

  // Force POST for mutations
  const finalMethod = isMutation ? "POST" : method

  let response: Response

  if (finalMethod === "GET") {
    // Normalize query to reduce URL length and ensure consistent caching
    const normalizedQuery = normalizeGraphQLQuery(query)

    // For GET, encode query and variables as URL parameters
    const params = new URLSearchParams()
    params.set("query", normalizedQuery)
    if (variables && Object.keys(variables).length > 0) {
      params.set("variables", JSON.stringify(variables))
    }

    response = await fetch(`/graphql/?${params.toString()}`, {
      method: "GET",
      credentials: "same-origin",
    })
  } else {
    // For POST, send in request body
    response = await fetch("/graphql/", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      credentials: "same-origin",
      body: JSON.stringify({
        query,
        variables,
      }),
    })
  }

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

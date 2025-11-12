import type { GraphQLResponse } from "../types"

function normalizeGraphQLQuery(query: string): string {
  return query
    .trim()
    .replace(/\s+/g, " ")
    .replace(/\s*([{}(),:])\s*/g, "$1")
}

export async function fetchGraphQL<T>(
  query: string,
  variables?: Record<string, unknown>,
  options?: { method?: "GET" | "POST"; operationName?: string }
): Promise<T> {
  const method = options?.method ?? "GET"
  const isMutation = query.trim().startsWith("mutation")
  const finalMethod = isMutation ? "POST" : method

  let response: Response

  if (finalMethod === "GET") {
    const params = new URLSearchParams()
    params.set("query", normalizeGraphQLQuery(query))
    if (variables && Object.keys(variables).length > 0) {
      params.set("variables", JSON.stringify(variables))
    }
    // Send operation name for ETag support
    if (options?.operationName) {
      params.set("operationName", options.operationName)
    }
    response = await fetch(`/graphql/?${params.toString()}`, {
      method: "GET",
      credentials: "same-origin",
    })
  } else {
    response = await fetch("/graphql/", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      credentials: "same-origin",
      body: JSON.stringify({ query, variables }),
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

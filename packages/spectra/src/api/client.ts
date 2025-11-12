import type { GraphQLResponse } from "../types"

function normalizeGraphQLQuery(query: string): string {
  return query
    .trim()
    .replace(/\s+/g, " ")
    .replace(/\s*([{}(),:])\s*/g, "$1")
}

/**
 * Detect if we're running in Safari/WebKit
 * Safari's fetch() API doesn't automatically send If-None-Match headers,
 * even though the browser cache exists. We only use manual ETag handling for Safari.
 */
function isSafari(): boolean {
  // Check for Safari-specific features
  // @ts-expect-error - Safari-specific API
  return (
    typeof window !== "undefined" &&
    window.safari !== undefined &&
    window.safari.pushNotification !== undefined
  )
}

/**
 * Manual ETag cache using localStorage (Safari only)
 *
 * Safari/WebKit fetch() doesn't automatically send If-None-Match headers,
 * even in Safari 17.2+. This is specific to the fetch() API, not browser navigation.
 *
 * For Safari only, we manually:
 * 1. Store ETags from responses in localStorage
 * 2. Include stored ETags as If-None-Match headers in subsequent requests
 * 3. Handle 304 responses by returning cached data from localStorage
 *
 * Other browsers use native HTTP cache behavior (respects hard refresh).
 */
const ETAG_CACHE_PREFIX = "fpbase_etag_"
const DATA_CACHE_PREFIX = "fpbase_data_"

function getStoredETag(url: string): string | null {
  try {
    return localStorage.getItem(ETAG_CACHE_PREFIX + url)
  } catch {
    return null
  }
}

function storeETag(url: string, etag: string): void {
  try {
    localStorage.setItem(ETAG_CACHE_PREFIX + url, etag)
  } catch {
    // Ignore storage errors
  }
}

function getStoredData<T>(url: string): T | null {
  try {
    const cached = localStorage.getItem(DATA_CACHE_PREFIX + url)
    if (!cached) {
      console.warn(`[ETag Cache] No stored data found for ${url}`)
      return null
    }
    return JSON.parse(cached)
  } catch (error) {
    console.error(`[ETag Cache] Failed to parse cached data:`, error)
    return null
  }
}

function storeData(url: string, data: unknown): void {
  try {
    const serialized = JSON.stringify(data)
    localStorage.setItem(DATA_CACHE_PREFIX + url, serialized)
    console.log(`[ETag Cache] Stored ${(serialized.length / 1024).toFixed(2)}KB for ${url}`)
  } catch (error) {
    console.error(`[ETag Cache] Failed to store data:`, error)
  }
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
  let url: string | undefined

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

    url = `/graphql/?${params.toString()}`

    // const useSafariWorkaround = isSafari()
    const useSafariWorkaround = false // for TESTING

    if (useSafariWorkaround) {
      // Safari: Manual ETag handling with localStorage
      const storedETag = getStoredETag(url)
      const headers: HeadersInit = {}
      if (storedETag) {
        headers["If-None-Match"] = storedETag
        console.log(`[Safari ETag] Sending If-None-Match: ${storedETag}`)
      }

      response = await fetch(url, {
        method: "GET",
        credentials: "same-origin",
        headers,
        cache: "no-store", // Disable browser cache for manual handling
      })

      // Handle 304 Not Modified
      if (response.status === 304) {
        console.log(`[Safari ETag] Got 304, retrieving from localStorage`)
        const cachedData = getStoredData<T>(url)
        if (cachedData) {
          console.log(`[Safari ETag] Returning cached data`)
          return cachedData
        }
        // Fallback: refetch without ETag
        console.warn(`[Safari ETag] Got 304 but no cached data - refetching`)
        response = await fetch(url, {
          method: "GET",
          credentials: "same-origin",
          cache: "no-store",
        })
      }

      // Store new ETag and data
      const etag = response.headers.get("etag")
      if (etag) {
        storeETag(url, etag)
      }
    } else {
      // Chrome/Firefox: Use native browser HTTP cache (respects hard refresh)
      response = await fetch(url, {
        method: "GET",
        credentials: "same-origin",
        // Use default cache mode - browser handles ETags automatically
      })
    }
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

  // Store data in localStorage for Safari ETag cache (GET requests only)
  if (finalMethod === "GET" && url && isSafari()) {
    storeData(url, json.data)
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

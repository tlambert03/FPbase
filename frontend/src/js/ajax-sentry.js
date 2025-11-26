/**
 * Fetch Error Tracking for Sentry
 *
 * This module provides centralized error tracking for fetch requests.
 * It automatically captures failed requests and reports them to Sentry with full context.
 *
 * Automatically tracks:
 * - Network failures
 * - HTTP error responses (4xx, 5xx)
 * - Timeout errors
 * - Parse errors
 *
 * Usage:
 *   Use fetchWithSentry() instead of fetch()
 *
 * Noise reduction hooks:
 *   - Set skipSentry: true in options to skip reporting
 *   - Customize shouldReportError() for global filtering
 */

import * as Sentry from "@sentry/browser"

/**
 * Determine if an error should be reported to Sentry
 * Override this function to customize filtering logic
 */
export function shouldReportError(_url, status, _method) {
  // Don't report user aborts
  if (status === 0) return false

  // Don't report 503 Service Unavailable
  // these are transient infrastructure issues that resolve on retry
  if (status === 503) return false

  // Add more filtering rules here as needed
  // Example: Don't report 404s for optional resources
  // if (status === 404 && url.includes('/optional/')) return false
  // Example: Filter by method
  // if (_method === 'GET' && status === 404) return false

  return true
}

/**
 * Fetch wrapper with automatic Sentry error reporting
 *
 * Drop-in replacement for fetch() that automatically reports errors to Sentry.
 *
 * Usage:
 *   import { fetchWithSentry } from './ajax-sentry'
 *   fetchWithSentry('/api/endpoint', { method: 'POST', ... })
 *
 * Options:
 *   - skipSentry: true - Skip Sentry reporting for this request
 *   - All standard fetch options are supported
 *
 * @param {string} url - The URL to fetch
 * @param {object} options - Fetch options (method, headers, body, etc.)
 * @returns {Promise<Response>} - The fetch response
 */
export async function fetchWithSentry(url, options = {}) {
  const method = options.method || "GET"
  const startTime = Date.now()

  try {
    const response = await fetch(url, options)

    // Check if response is an error status
    if (!response.ok) {
      // Apply filtering logic
      if (shouldReportError(url, response.status, method) && !options.skipSentry) {
        const duration = Date.now() - startTime
        const errorMessage = `${response.status} ${response.statusText}`

        // Try to get response body for context (clone so original can still be used)
        let responseText
        try {
          const cloned = response.clone()
          responseText = await cloned.text()
        } catch {
          // Ignore errors when trying to read response body
          responseText = undefined
        }

        // Create detailed error for Sentry
        const error = new Error(`Fetch ${method} ${url} failed: ${errorMessage}`)
        error.name = "FetchError"

        // Capture to Sentry with full context
        Sentry.captureException(error, {
          level: response.status >= 500 ? "error" : "warning",
          tags: {
            ajax_url: url,
            ajax_method: method,
            http_status: response.status,
            error_type: "fetch",
          },
          contexts: {
            fetch: {
              url,
              method,
              status: response.status,
              statusText: response.statusText,
              responseText: responseText ? responseText.substring(0, 1000) : undefined,
              duration,
              headers: Object.fromEntries(response.headers.entries()),
              requestHeaders: options.headers || {},
            },
          },
          fingerprint: ["fetch", method, url, String(response.status)],
        })

        // Log to console for debugging
        if (process.env.NODE_ENV !== "production") {
          console.error("Fetch Error:", {
            url,
            method,
            status: response.status,
            statusText: response.statusText,
            duration,
            response: responseText,
          })
        }
      }
    }

    return response
  } catch (error) {
    // Network error, timeout, or other fetch failure
    const duration = Date.now() - startTime

    // Apply filtering logic (status 0 for network errors)
    if (shouldReportError(url, 0, method) && !options.skipSentry) {
      // Create detailed error for Sentry
      const sentryError = new Error(`Fetch ${method} ${url} network error: ${error.message}`)
      sentryError.name = "FetchNetworkError"

      // Capture to Sentry with full context
      Sentry.captureException(sentryError, {
        level: "error",
        tags: {
          ajax_url: url,
          ajax_method: method,
          http_status: 0,
          error_type: "fetch_network",
        },
        contexts: {
          fetch: {
            url,
            method,
            error: error.message,
            duration,
            requestHeaders: options.headers || {},
          },
        },
        fingerprint: ["fetch", method, url, "network_error"],
      })

      // Log to console for debugging
      if (process.env.NODE_ENV !== "production") {
        console.error("Fetch Network Error:", {
          url,
          method,
          error: error.message,
          duration,
        })
      }
    }

    // Re-throw the error so calling code can handle it
    throw error
  }
}

// Export fetchWithSentry to global window object for use in inline scripts
if (typeof window !== "undefined") {
  window.fetchWithSentry = fetchWithSentry
}

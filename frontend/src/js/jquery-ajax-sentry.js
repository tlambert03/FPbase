/**
 * Global jQuery AJAX Error Tracking for Sentry
 *
 * This module sets up global AJAX error handlers to automatically capture
 * failed requests and report them to Sentry with full context.
 *
 * Automatically tracks:
 * - Network failures
 * - HTTP error responses (4xx, 5xx)
 * - Timeout errors
 * - Parse errors
 *
 * Import this in any entry point that uses jQuery AJAX calls.
 */

import * as Sentry from "@sentry/browser"
import $ from "jquery"

/**
 * Setup global jQuery AJAX error tracking
 */
export function setupAjaxErrorTracking() {
  // Track all AJAX errors globally
  $(document).ajaxError((_event, jqXHR, ajaxSettings, thrownError) => {
    // Don't report if it's a user abort
    if (jqXHR.statusText === "abort") {
      return
    }

    // Construct error message
    const errorMessage = thrownError || jqXHR.statusText || "Unknown AJAX Error"
    const url = ajaxSettings.url || "unknown URL"
    const method = ajaxSettings.type || ajaxSettings.method || "GET"

    // Create detailed error for Sentry
    const error = new Error(`AJAX ${method} ${url} failed: ${errorMessage}`)
    error.name = "AjaxError"

    // Capture to Sentry with full context
    Sentry.captureException(error, {
      level: jqXHR.status >= 500 ? "error" : "warning",
      tags: {
        ajax_url: url,
        ajax_method: method,
        http_status: jqXHR.status || "unknown",
        error_type: "ajax",
      },
      contexts: {
        ajax: {
          url,
          method,
          status: jqXHR.status,
          statusText: jqXHR.statusText,
          responseText: jqXHR.responseText ? jqXHR.responseText.substring(0, 1000) : undefined,
          readyState: jqXHR.readyState,
          settings: {
            contentType: ajaxSettings.contentType,
            dataType: ajaxSettings.dataType,
            async: ajaxSettings.async,
            timeout: ajaxSettings.timeout,
          },
        },
      },
      fingerprint: ["ajax", method, url, String(jqXHR.status)],
    })

    // Log to console for debugging
    if (process.env.NODE_ENV !== "production") {
      console.error("AJAX Error:", {
        url,
        method,
        status: jqXHR.status,
        error: errorMessage,
        response: jqXHR.responseText,
      })
    }
  })
}

// Auto-setup on import
setupAjaxErrorTracking()

export default setupAjaxErrorTracking

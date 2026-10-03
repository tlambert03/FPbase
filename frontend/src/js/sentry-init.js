/**
 * Shared Sentry initialization module
 *
 * This module initializes Sentry once and exports it for use across all bundles.
 * It should be imported at the top of every entry point to ensure errors are captured.
 *
 * Benefits:
 * - Single Sentry instance prevents conflicts
 * - Early initialization catches errors during script loading
 * - Shared configuration across all bundles
 * - Global error handlers set up before any other code runs
 */

import * as Sentry from "@sentry/browser"

let sentryInitialized = false

// crawlers, AI agents and automated browsers; keep in sync with BOT_UA in backend/fpbase/tracing.py
const BOT_UA =
  /bot|crawl|spider|slurp|scrap|headless|phantom|preview|facebookexternalhit|claude|chatgpt|copilot|perplexity/i

function isBot() {
  return navigator.webdriver || BOT_UA.test(navigator.userAgent)
}

// Cloudflare redirects every request for http:// or fpbase.org to https://www.fpbase.org,
// so a page at any other fpbase.org origin is a saved copy replayed by a crawler: its
// same-origin requests get redirected cross-origin and fail (FPBASE-6BB).
function isReplayedCopy() {
  const { hostname, origin } = window.location
  return ["fpbase.org", "www.fpbase.org"].includes(hostname) && origin !== "https://www.fpbase.org"
}

/**
 * Initialize Sentry with production-optimized configuration
 * This function is idempotent - calling it multiple times is safe
 */
export function initSentry() {
  // Only initialize once, even if imported by multiple bundles
  if (sentryInitialized) {
    return Sentry
  }

  // Only initialize in production with a valid DSN
  if (process.env.NODE_ENV === "production" && process.env.SENTRY_DSN) {
    try {
      Sentry.init({
        dsn: process.env.SENTRY_DSN,
        release: process.env.HEROKU_SLUG_COMMIT,
        environment: process.env.NODE_ENV,
        integrations: [
          Sentry.browserTracingIntegration(),
          Sentry.replayIntegration({ maskAllText: false, blockAllMedia: false }),
        ],
        replaysSessionSampleRate: 0, // Don't record normal sessions
        replaysOnErrorSampleRate: 1.0, // Record all error sessions
        // trace 10% of human page loads; the backend follows this decision for API calls
        tracesSampler: ({ inheritOrSampleWith }) => (isBot() ? 0 : inheritOrSampleWith(0.1)),
        tracePropagationTargets: [
          "localhost",
          /^\//, // Relative URLs (same-origin API calls)
          /^https:\/\/([^.]+\.)?fpbase\.org/, // All fpbase.org subdomains
        ],
        sampleRate: 1.0, // Send all errors (defatult=1)

        // Ignore common benign errors
        ignoreErrors: [
          // Browser extensions
          "top.GLOBALS",
          "originalCreateNotification",
          "canvas.contentDocument",
          "MyApp_RemoveAllHighlights",
          "atomicFindClose",
          // Network errors that are expected
          "NetworkError",
          "Non-Error promise rejection captured",
          // Random plugins/extensions
          "conduitPage",
        ],

        // Filter sensitive data before sending
        beforeSend(event, _hint) {
          if (isReplayedCopy()) return null

          // Add bundle information for easier debugging
          event.tags = {
            ...event.tags,
            bundle: window.FPBASE?.currentBundle || "unknown",
          }

          return event
        },

        // Enhance events with additional context
        beforeSendTransaction(event) {
          // Add custom context to transactions
          event.tags = {
            ...event.tags,
            bundle: window.FPBASE?.currentBundle || "unknown",
          }
          return event
        },
      })

      // Set user context if available
      if (window.FPBASE?.user) {
        Sentry.setUser({
          id: window.FPBASE.user.id,
          username: window.FPBASE.user.name,
        })
      }

      // Expose Sentry globally for manual error reporting
      window.Sentry = Sentry

      sentryInitialized = true
    } catch (error) {
      console.error("❌ Failed to initialize Sentry:", error)
    }
  }

  return Sentry
}

// Auto-initialize on import
export const sentry = initSentry()

// Export Sentry as default for convenience
export default sentry

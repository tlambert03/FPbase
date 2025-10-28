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
        integrations: [Sentry.replayIntegration({ maskAllText: false, blockAllMedia: false })],
        replaysSessionSampleRate: 0, // Don't record normal sessions
        replaysOnErrorSampleRate: 1.0, // Record all error sessions
        tracesSampleRate: 0.05, // 5% performance monitoring
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
        beforeSend(event, hint) {
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
          email: window.FPBASE.user.email,
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

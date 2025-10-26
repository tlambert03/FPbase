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

                // Session Replay - capture user sessions for debugging
                integrations: [
                    Sentry.replayIntegration({
                        maskAllText: false, // Capture text for better debugging
                        blockAllMedia: false, // Capture media elements
                    }),
                    // Browser profiling for performance insights
                    Sentry.browserProfilingIntegration(),
                ],

                // Session Replay sampling
                replaysSessionSampleRate: 0.1, // 10% of normal sessions
                replaysOnErrorSampleRate: 1.0, // 100% of error sessions

                // Performance monitoring
                tracesSampleRate: 0.1, // 10% of transactions
                profilesSampleRate: 0.1, // 10% of transactions for profiling

                // Trace propagation to backend APIs
                tracePropagationTargets: [
                    "fpbase.org",
                    /^https:\/\/.*\.fpbase\.org/,
                    /^\/api\//,
                    /^\/graphql\//,
                    /^\/ajax\//,
                ],

                // Privacy settings
                sendDefaultPii: false, // Don't send PII by default

                // Error sampling - capture all errors
                sampleRate: 1.0,

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
                    // Remove sensitive query parameters
                    if (event.request?.url) {
                        const url = new URL(event.request.url, window.location.origin)
                        const sensitiveParams = ["token", "key", "password", "secret", "api_key"]
                        sensitiveParams.forEach((param) => {
                            if (url.searchParams.has(param)) {
                                url.searchParams.set(param, "[Filtered]")
                            }
                        })
                        event.request.url = url.toString()
                    }

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

            console.log("✅ Sentry initialized successfully")
        } catch (error) {
            console.error("❌ Failed to initialize Sentry:", error)
        }
    } else if (process.env.NODE_ENV !== "production") {
        console.log("ℹ️ Sentry not initialized (development mode)")
    } else {
        console.warn("⚠️ Sentry not initialized (missing DSN)")
    }

    return Sentry
}

// Auto-initialize on import
export const sentry = initSentry()

// Export Sentry as default for convenience
export default sentry

/**
 * React Error Boundary with Sentry integration
 *
 * Catches errors in React component trees and reports them to Sentry.
 * Displays a fallback UI when errors occur instead of crashing the entire app.
 *
 * Usage:
 *   import { SentryErrorBoundary } from './js/sentry-error-boundary'
 *
 *   <SentryErrorBoundary fallback={<ErrorFallback />}>
 *     <YourApp />
 *   </SentryErrorBoundary>
 */

import { Component } from "react"
import * as Sentry from "@sentry/browser"

/**
 * Default fallback component shown when an error occurs
 */
function DefaultErrorFallback({ error, errorInfo, resetError }) {
    return (
        <div
            style={{
                padding: "2rem",
                margin: "2rem auto",
                maxWidth: "600px",
                border: "2px solid #dc3545",
                borderRadius: "8px",
                backgroundColor: "#f8d7da",
                color: "#721c24",
            }}
        >
            <h2 style={{ marginTop: 0 }}>
                <i className="fas fa-exclamation-triangle" style={{ marginRight: "0.5rem" }}></i>
                Something Went Wrong
            </h2>
            <p>We're sorry, but something unexpected happened. Our team has been notified.</p>

            {error && (
                <details style={{ marginTop: "1rem", fontSize: "0.9em" }}>
                    <summary style={{ cursor: "pointer", fontWeight: "bold" }}>
                        Error Details
                    </summary>
                    <pre
                        style={{
                            marginTop: "0.5rem",
                            padding: "1rem",
                            backgroundColor: "#fff",
                            border: "1px solid #ccc",
                            borderRadius: "4px",
                            overflow: "auto",
                            fontSize: "0.85em",
                        }}
                    >
                        {error.toString()}
                        {errorInfo?.componentStack}
                    </pre>
                </details>
            )}

            {resetError && (
                <button
                    onClick={resetError}
                    style={{
                        marginTop: "1rem",
                        padding: "0.5rem 1rem",
                        backgroundColor: "#007bff",
                        color: "#fff",
                        border: "none",
                        borderRadius: "4px",
                        cursor: "pointer",
                        fontSize: "1rem",
                    }}
                >
                    Try Again
                </button>
            )}
        </div>
    )
}

/**
 * Error Boundary component that catches React errors and reports to Sentry
 */
export class SentryErrorBoundary extends Component {
    constructor(props) {
        super(props)
        this.state = {
            hasError: false,
            error: null,
            errorInfo: null,
            eventId: null,
        }
    }

    static getDerivedStateFromError(error) {
        // Update state so the next render will show the fallback UI
        return { hasError: true, error }
    }

    componentDidCatch(error, errorInfo) {
        // Log error to Sentry
        const eventId = Sentry.captureException(error, {
            contexts: {
                react: {
                    componentStack: errorInfo.componentStack,
                },
            },
            tags: {
                errorBoundary: this.props.name || "SentryErrorBoundary",
            },
        })

        this.setState({
            error,
            errorInfo,
            eventId,
        })

        // Also log to console in development
        if (process.env.NODE_ENV !== "production") {
            console.error("React Error Boundary caught an error:", error, errorInfo)
        }
    }

    resetError = () => {
        this.setState({
            hasError: false,
            error: null,
            errorInfo: null,
            eventId: null,
        })

        if (this.props.onReset) {
            this.props.onReset()
        }
    }

    render() {
        if (this.state.hasError) {
            // Use custom fallback if provided, otherwise use default
            const FallbackComponent = this.props.fallback || DefaultErrorFallback

            if (typeof FallbackComponent === "function") {
                return (
                    <FallbackComponent
                        error={this.state.error}
                        errorInfo={this.state.errorInfo}
                        eventId={this.state.eventId}
                        resetError={this.resetError}
                    />
                )
            }

            // If fallback is a React element, render it directly
            return FallbackComponent
        }

        return this.props.children
    }
}

/**
 * Convenience function to wrap a component with error boundary
 *
 * Usage:
 *   export default withSentryErrorBoundary(MyComponent, { name: 'MyComponent' })
 */
export function withSentryErrorBoundary(Component, options = {}) {
    const WrappedComponent = (props) => (
        <SentryErrorBoundary {...options}>
            <Component {...props} />
        </SentryErrorBoundary>
    )

    WrappedComponent.displayName = `withSentryErrorBoundary(${Component.displayName || Component.name || "Component"})`

    return WrappedComponent
}

export default SentryErrorBoundary

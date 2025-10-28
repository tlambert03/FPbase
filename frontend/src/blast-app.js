// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"

import App from "@fpbase/blast"
import React from "react"
import { createRoot } from "react-dom/client"
import { SentryErrorBoundary } from "./js/sentry-error-boundary"

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE.currentBundle = "blast"

const root = createRoot(document.getElementById("blast-app"))
root.render(
  <SentryErrorBoundary name="BlastApp">
    <App />
  </SentryErrorBoundary>
)

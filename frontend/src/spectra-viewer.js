// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"

import "./js/detect-touch"
import { createElement } from "react"
import { createRoot } from "react-dom/client"
import App from "@fpbase/spectra"
import { SentryErrorBoundary } from "./js/sentry-error-boundary"

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE.currentBundle = "spectraViewer"

const root = createRoot(document.getElementById("spectra-viewer"))
root.render(
    createElement(SentryErrorBoundary, { name: "SpectraViewer" },
        createElement(App, { uri: "/graphql/" }, null)
    )
)

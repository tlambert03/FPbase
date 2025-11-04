// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"

import "./js/detect-touch"
import { SimpleSpectraViewer } from "@fpbase/spectra"
import { createElement } from "react"
import { createRoot } from "react-dom/client"
import { SentryErrorBoundary } from "./js/sentry-error-boundary"

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE.currentBundle = "simpleSpectraViewer"

const elem = document.getElementById("spectra-viewer")

window.addEventListener("load", () => {
  const root = createRoot(elem)

  // Check if data attributes are present (protein_detail.html)
  const hasDataAttrs = elem.hasAttribute("data-spectra")

  let props
  if (hasDataAttrs) {
    // Use data attributes from template (protein_detail.html, compare.html)
    // CRITICAL: This ensures protein pages always show their specific spectra,
    // ignoring any URL params or session storage

    // Parse and ensure IDs are strings (Django templates may output numbers)
    const activeSpectra = JSON.parse(elem.getAttribute("data-spectra")).map(String)
    const hiddenSpectra = (JSON.parse(elem.getAttribute("data-hidden")) || []).map(String)

    props = {
      state: {
        activeSpectra,
        hiddenSpectra,
        chartOptions: JSON.parse(elem.getAttribute("data-options")),
      },
      fromUrl: false,
    }
  } else {
    // No data attributes - load from URL (spectra_graph.html)
    props = {
      fromUrl: true,
    }
  }

  root.render(
    createElement(
      SentryErrorBoundary,
      { name: "SimpleSpectraViewer" },
      createElement(SimpleSpectraViewer, props)
    )
  )
})

// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"

import mount from "@fpbase/protein-table"

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE.currentBundle = "proteinTable"

// Mount the protein table app when the DOM is ready
if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", () => {
    const container = document.getElementById("protein-table")
    if (container) {
      try {
        mount(container)
      } catch (error) {
        console.error("Failed to mount protein table:", error)
        if (window.Sentry) {
          window.Sentry.captureException(error)
        }
      }
    }
  })
} else {
  const container = document.getElementById("protein-table")
  if (container) {
    try {
      mount(container)
    } catch (error) {
      console.error("Failed to mount protein table:", error)
      if (window.Sentry) {
        window.Sentry.captureException(error)
      }
    }
  }
}

// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/ajax-sentry.js" // Track AJAX and fetch errors

// FRET calculator functionality with Highcharts
import initFRET from "./js/fret.js"

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE = {
  ...window.FPBASE,
  currentBundle: "fret",
  initFRET,
}

// Also expose initFRET globally for legacy inline scripts
window.initFRET = initFRET

// Auto-initialize FRET calculator when bundle loads
document.addEventListener("DOMContentLoaded", () => {
  try {
    initFRET()
  } catch (error) {
    console.error("Error initializing FRET:", error)
    if (window.Sentry) {
      window.Sentry.captureException(error, {
        tags: { component: "fret-init" },
      })
    }
  }
})

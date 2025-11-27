/**
 * Enhanced Spectrum Submission Form V2
 *
 * Entry point for the new spectrum submission form with:
 * - Multi-column CSV/TSV file support
 * - Column picker UI
 * - Client-side normalization
 * - Interactive Highcharts preview with range selector
 */

// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/ajax-sentry.js"
import { initSpectrumForm } from "./js/spectrum-form/form-controller.js"

// Initialize when DOM is ready
document.addEventListener("DOMContentLoaded", initSpectrumForm)

// Expose for debugging
window.FPBASE = window.FPBASE || {}
window.FPBASE.spectrumFormV2 = {
  reinit: initSpectrumForm,
}

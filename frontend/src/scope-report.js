// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/ajax-sentry.js" // Track jQuery AJAX errors

// Scope report functionality with Highcharts
import "./js/scope_report.js"

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE.currentBundle = "scope-report"

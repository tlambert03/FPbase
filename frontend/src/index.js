// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/jquery-ajax-sentry.js" // Track jQuery AJAX errors

import "regenerator-runtime/runtime"

// jQuery, select2, and nouislider loaded from CDN in base.html
import "select2/dist/css/select2.css"
import "select2-theme-bootstrap4/dist/select2-bootstrap.css"
import "nouislider/distribute/nouislider.min.css"
import "./css/style.scss"

// Bootstrap loaded from CDN in base.html for jQuery compatibility

import "./js/project.js"
import initSearch from "./js/search_logic.js"
import "./js/protein_page.js"
import "./js/favit.js"
import "./js/jquery.formset.js"
import "./js/onload.js"
// microscope.js loaded separately via CDN on microscope pages
import "./js/scope_report.js"

import * as d3 from "d3"
import initAutocomplete from "./js/algolia.js"
import initFRET from "./js/fret.js"
import FPPropChart from "./js/ichart.js"
import LineageChart from "./js/lineage.js"

// Mark this bundle for Sentry context
// Save D3 v7 to a separate variable and set as default if no other d3 is loaded
// This allows microscope pages to load D3 v3 from CDN without conflict
window.d3v7 = d3
if (!window.d3) {
  window.d3 = d3
}

window.FPBASE = window.FPBASE || {}
window.FPBASE = {
  ...window.FPBASE,
  currentBundle: "main",
  initAutocomplete,
  initSearch,
  FPPropChart,
  LineageChart,
  initFRET,
}

// Also expose initSearch globally for legacy inline scripts
window.initSearch = initSearch

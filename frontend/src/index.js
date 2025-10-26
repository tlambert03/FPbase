// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/jquery-ajax-sentry.js" // Track jQuery AJAX errors

import "regenerator-runtime/runtime"
import "select2/dist/css/select2.css"
import "select2-theme-bootstrap4/dist/select2-bootstrap.css"
import "nouislider/distribute/nouislider.min.css"
import "./css/style.scss"

import "bootstrap"

import "select2/dist/js/select2.full.js"

import "./js/project.js"
import "./js/search_logic.js"
import "./js/protein_page.js"
import "./js/favit.js"
import "./js/jquery.formset.js"
import "./js/onload.js"
// microscope.js loaded separately via CDN on microscope pages
import "./js/scope_report.js"

import FPPropChart from "./js/ichart.js"
import initAutocomplete from "./js/algolia.js"
import LineageChart from "./js/lineage.js"
import initFRET from "./js/fret.js"
import * as d3 from "d3";

// Mark this bundle for Sentry context
window.d3 = d3;

window.FPBASE = window.FPBASE || {}
window.FPBASE = {
  ...window.FPBASE,
  currentBundle: "main",
  initAutocomplete,
  FPPropChart,
  LineageChart,
  initFRET
}

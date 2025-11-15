// NOTE:
// jQuery, select2, and nouislider loaded from CDN in base.html
// Bootstrap loaded from CDN in base.html for jQuery compatibility
// microscope.js loaded separately via CDN on microscope pages
// scope_report.js and fret.js moved to separate bundles to avoid loading Highcharts on all pages
// D3 and chart components moved to separate bundle (d3-charts.js) to avoid loading D3 on all pages

// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/jquery-ajax-sentry.js" // Track jQuery AJAX errors
import { icon } from "./js/icons.js" // Icon helper for dynamic HTML

import "select2/dist/css/select2.css"
import "select2-theme-bootstrap4/dist/select2-bootstrap.css"
// NOTE: nouislider CSS moved to d3-charts.js since it's only used on ichart page
import "./css/style.scss"

import "./js/project.js"
import initSearch from "./js/search_logic.js"
import "./js/favit.js"
import "./js/jquery.formset.js"
import "./js/onload.js"

// D3 and chart components removed - now lazy-loaded via d3-charts.js
import initAutocomplete from "./js/algolia.js"

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE = {
  ...window.FPBASE,
  currentBundle: "main",
  initAutocomplete,
  initSearch,

  // Lazy-load D3 charts only when needed
  _d3ChartsPromise: null,
  async loadD3Charts() {
    if (!this._d3ChartsPromise) {
      this._d3ChartsPromise = import("./d3-charts.js")
        .then((module) => {
          this.FPPropChart = module.FPPropChart
          this.LineageChart = module.LineageChart
          return module
        })
        .catch((error) => {
          console.error("Failed to load D3 charts:", error)
          if (window.Sentry) {
            window.Sentry.captureException(error, {
              tags: { component: "d3-charts-lazy-load" },
            })
          }
          throw error
        })
    }
    return this._d3ChartsPromise
  },
}

// Also expose initSearch globally for legacy inline scripts
window.initSearch = initSearch

// Notify that FPBASE is ready for use by inline scripts
// This is needed because vite bundles load with defer, so inline scripts
// at the end of <body> may execute before this module
window.dispatchEvent(
  new CustomEvent("fpbase:ready", {
    bubbles: true,
    detail: { FPBASE: window.FPBASE },
  })
)

// Initialize autocomplete search when DOM is ready
document.addEventListener("DOMContentLoaded", () => {
  // Initialize autocomplete if search input exists
  if (document.getElementById("algolia-search-input")) {
    initAutocomplete()
  }

  // Auto-initialization: Look for elements with data-fpbase-init attribute
  document.querySelectorAll("[data-fpbase-init]").forEach((element) => {
    const initType = element.dataset.fpbaseInit

    try {
      switch (initType) {
        case "search": {
          const fields = JSON.parse(element.dataset.filterFields || "{}")
          const operators = JSON.parse(element.dataset.filterOperators || "{}")
          const labels = JSON.parse(element.dataset.filterLabels || "{}")
          initSearch(fields, operators, labels)
          break
        }
        case "molstar": {
          // Lazy load Mol* only when the structure section becomes visible
          const pdbIds = JSON.parse(element.dataset.pdbIds || "[]")

          if (!pdbIds || pdbIds.length === 0) {
            console.warn("Mol* init: No PDB IDs provided")
            break
          }

          // Use Intersection Observer for lazy loading
          const observer = new IntersectionObserver(
            (entries) => {
              entries.forEach((entry) => {
                if (entry.isIntersecting) {
                  // Dynamically import Mol* bundle only when needed
                  import("./my-molstar.js")
                    .then((module) => {
                      const initPDB = module.default
                      initPDB(pdbIds)
                    })
                    .catch((error) => {
                      console.error("Failed to load Mol*:", error)
                      if (window.Sentry) {
                        window.Sentry.captureException(error, {
                          tags: { component: "molstar-lazy-load" },
                        })
                      }

                      // Show user-friendly error message
                      const links = pdbIds
                        .map((id) => `<a href="https://www.rcsb.org/structure/${id}">${id}</a>`)
                        .join(", ")

                      element.innerHTML = `
                        <div class="col-12">
                          <div class="alert alert-warning" role="alert">
                            <h5 class="alert-heading">
                              ${icon("alert", "mr-2")}
                              Unable to Load 3D Structure Viewer
                            </h5>
                            <p class="mb-2">The molecular structure viewer failed to load. This may be due to a network issue or browser compatibility problem.</p>
                            <hr>
                            <p class="mb-0">
                              You can view these structures directly at RCSB PDB: ${links}
                            </p>
                          </div>
                        </div>
                      `
                    })

                  // Stop observing once loaded
                  observer.unobserve(entry.target)
                }
              })
            },
            {
              // Load when element is 200px from entering viewport
              rootMargin: "200px",
            }
          )

          observer.observe(element)
          break
        }
        case "microscope": {
          // Microscope.js is loaded as a legacy script (not ES module)
          // It defines window.initMicroscope when ready
          // Wait for both DOM and dependencies (d3, nvd3, jQuery) from CDN

          const initWhenReady = () => {
            if (
              typeof window.initMicroscope === "function" &&
              typeof nv !== "undefined" &&
              typeof d3 !== "undefined" &&
              typeof $ !== "undefined"
            ) {
              // Set global variables that microscope.js expects
              if (element.dataset.scopeSpectra) {
                window.scopespectra = JSON.parse(element.dataset.scopeSpectra)
              }
              if (element.dataset.scopeConfig) {
                window.scopecfg = JSON.parse(element.dataset.scopeConfig)
              }

              window.initMicroscope()
            } else {
              // Dependencies not ready yet, retry
              console.warn("Microscope dependencies not loaded yet, retrying...")
              setTimeout(initWhenReady, 500)
            }
          }

          initWhenReady()
          break
        }
        case "ichart": {
          // Lazy-load D3 charts and initialize interactive chart
          const apiUrl = element.dataset.apiUrl

          if (!apiUrl) {
            console.error("ichart init: No API URL provided")
            break
          }

          window.FPBASE.loadD3Charts()
            .then(({ FPPropChart }) => {
              return fetch(apiUrl, {
                headers: { Accept: "application/json" },
              })
                .then((response) => response.json())
                .then((data) => {
                  if (data == null || !Array.isArray(data)) {
                    throw new Error("Invalid data format from API")
                  }

                  // Initialize chart with data
                  const fpchart = FPPropChart().data(data)
                  window.d3v7.select("#graph").call(fpchart)

                  // Hide loading indicator
                  const loadingGif = element.querySelector(".loadinggif")
                  if (loadingGif) {
                    loadingGif.style.display = "none"
                  }
                })
            })
            .catch((error) => {
              console.error("Error loading interactive chart:", error)

              // Hide loading indicator
              const loadingGif = element.querySelector(".loadinggif")
              if (loadingGif) {
                loadingGif.style.display = "none"
              }

              // Show error modal if available
              const errorModal = document.getElementById("errorModal")
              if (errorModal && window.$ && window.$.fn.modal) {
                window.$(errorModal).modal("show")
              }

              if (window.Sentry) {
                window.Sentry.captureException(error, {
                  tags: { component: "ichart-init" },
                })
              }
            })
          break
        }
        default:
          console.warn(`Unknown init type: ${initType}`)
      }
    } catch (error) {
      console.error(`Error initializing ${initType}:`, error)
    }
  })
})

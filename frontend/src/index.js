// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/jquery-ajax-sentry.js" // Track jQuery AJAX errors

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

// Auto-initialization: Look for elements with data-fpbase-init attribute
document.addEventListener("DOMContentLoaded", () => {
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
        case "autocomplete":
          initAutocomplete()
          break
        case "fret":
          initFRET()
          break
        case "litemol": {
          // Lazy load LiteMol only when the structure section becomes visible
          const pdbIds = JSON.parse(element.dataset.pdbIds || "[]")

          if (!pdbIds || pdbIds.length === 0) {
            console.warn("LiteMol init: No PDB IDs provided")
            break
          }

          // Use Intersection Observer for lazy loading
          const observer = new IntersectionObserver(
            (entries) => {
              entries.forEach((entry) => {
                if (entry.isIntersecting) {
                  // Dynamically import LiteMol bundle only when needed
                  import("./my-litemol.js")
                    .then((module) => {
                      const initPDB = module.default
                      initPDB(pdbIds)
                    })
                    .catch((error) => {
                      console.error("Failed to load LiteMol:", error)
                      if (window.Sentry) {
                        window.Sentry.captureException(error, {
                          tags: { component: "litemol-lazy-load" },
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
                              <i class="fas fa-exclamation-triangle mr-2"></i>
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
        default:
          console.warn(`Unknown init type: ${initType}`)
      }
    } catch (error) {
      console.error(`Error initializing ${initType}:`, error)
    }
  })
})

// Dispatch custom event to signal bundle is ready
// This runs immediately - no need to wait for DOMContentLoaded
window.dispatchEvent(
  new CustomEvent("fpbase:ready", {
    detail: { bundle: "main", version: window.FPBASE },
  })
)

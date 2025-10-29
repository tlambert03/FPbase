// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"

import $ from "jquery"
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
  root.render(
    createElement(
      SentryErrorBoundary,
      { name: "SimpleSpectraViewer" },
      createElement(SimpleSpectraViewer, {
        ids: JSON.parse(elem.getAttribute("data-spectra")),
        options: JSON.parse(elem.getAttribute("data-options")),
        hidden: JSON.parse(elem.getAttribute("data-hidden")) || [],
      })
    )
  )
})
const name = elem.getAttribute("data-name")

if (name) {
  const svgElem = document.getElementById("simple_spectrasvg")
  $(svgElem).prepend(
    $("<desc>", { id: "svgDesc" }).text(
      `Fluorescent protein ${name} excitation and emission spectra`
    )
  )
  $(svgElem).prepend($("<title>", { id: "svgTitle" }).text(`${name} Spectrum`))
}

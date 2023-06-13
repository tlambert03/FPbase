import $ from "jquery"
import "./js/detect-touch"
import { createElement } from "react"
import { render } from "react-dom"
import { SimpleSpectraViewer } from "@fpbase/spectra"

const elem = document.getElementById("spectra-viewer")

window.onload = function() {
  render(
    createElement(SimpleSpectraViewer, {
      ids: JSON.parse(elem.getAttribute("data-spectra")),
      options: JSON.parse(elem.getAttribute("data-options")),
      hidden: JSON.parse(elem.getAttribute("data-hidden")) || [],
    }),
    elem
  )
}
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

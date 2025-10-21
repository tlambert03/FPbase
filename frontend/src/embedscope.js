import "jquery"
import "bootstrap"
import "nouislider"
import * as d3 from "d3"
import "select2/dist/js/select2.full.js"
// microscope.js loaded separately via CDN on microscope pages

// allow parent to apply styles to iframe
window.addEventListener("message", (event) => {
  if (event.data.type === "apply-css") {
    const style = document.createElement("style")
    style.innerHTML = event.data.css
    document.head.appendChild(style)
    console.log("iframe styles applied", style)
  }
})

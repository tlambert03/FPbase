// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"

import "jquery"
import "bootstrap"
import "nouislider"
// D3 v3 is loaded from CDN for microscope.js compatibility (see microscope_embed.html)
// Do NOT import D3 v7 here as it will overwrite the global d3 object and break nvd3
import "select2/dist/js/select2.full.js"
// microscope.js loaded separately via CDN on microscope pages

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE.currentBundle = "embedscope"

// allow parent to apply styles to iframe
window.addEventListener("message", (event) => {
  if (event.data.type === "apply-css") {
    const style = document.createElement("style")
    style.innerHTML = event.data.css
    document.head.appendChild(style)
    console.log("iframe styles applied", style)
  }
})

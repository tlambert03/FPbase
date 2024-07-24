import "jquery"
import "bootstrap"
import "nouislider"
import d3 from "d3"
import "./js/nv.d3.js"
import "select2/dist/js/select2.full.js"
import "./js/my-fontawesome.js"
import "./js/microscope.js"

// allow parent to apply styles to iframe
window.addEventListener("message", (event) => {
  if (event.data.type === "apply-css") {
    const style = document.createElement("style")
    style.innerHTML = event.data.css
    document.head.appendChild(style)
    console.log("iframe styles applied", style)
  }
})

/**
 * D3 Charts Bundle
 *
 * This bundle contains D3.js and chart components that are only used on specific pages.
 * It's lazy-loaded to avoid bloating the main bundle.
 *
 * Used on: ichart.html, lineage.html, protein_detail.html, organism_detail.html
 */

import * as d3 from "d3"
import noUiSlider from "nouislider"
import "nouislider/distribute/nouislider.min.css"

import FPPropChart from "./js/ichart.js"
import LineageChart from "./js/lineage.js"

// Save D3 v7 to a separate variable and set as default if no other d3 is loaded
// This allows microscope pages to load D3 v3 from CDN without conflict
window.d3v7 = d3
if (!window.d3) {
  window.d3 = d3
}

// Export noUiSlider for ichart.js which uses it
export { FPPropChart, LineageChart, noUiSlider }
export default { FPPropChart, LineageChart, noUiSlider }

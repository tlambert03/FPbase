/**
 * D3 Charts Bundle
 *
 * This bundle contains D3.js and chart components that are only used on specific pages.
 * It's lazy-loaded to avoid bloating the main bundle.
 *
 * Used on: ichart.html, lineage.html, protein_detail.html, organism_detail.html
 *
 * NOTE: Uses modular D3 imports instead of `import * as d3` for better tree-shaking.
 * This reduces the bundle size by importing only the D3 modules we actually use.
 */

import noUiSlider from "nouislider"
import "nouislider/distribute/nouislider.min.css"

import { max, min } from "d3-array"
import { axisBottom, axisLeft, axisRight, axisTop } from "d3-axis"
import { hsl } from "d3-color"
import { drag } from "d3-drag"
import { cluster, hierarchy, tree } from "d3-hierarchy"
import { scaleLinear, scaleLog } from "d3-scale"
// Re-export the D3 modules used by the charts for window.d3 (for template code compatibility)
import { select, selectAll } from "d3-selection"
import { zoom } from "d3-zoom"
// Import chart modules - they use modular D3 imports for tree-shaking
import FPPropChart from "./js/ichart.js"
import LineageChart from "./js/lineage.js"

// Create a minimal d3 namespace for backward compatibility with template code
const d3 = {
  select,
  selectAll,
  scaleLinear,
  scaleLog,
  axisBottom,
  axisLeft,
  axisTop,
  axisRight,
  hierarchy,
  tree,
  cluster,
  hsl,
  max,
  min,
  zoom,
  drag,
}

// Save D3 v7 to a separate variable and set as default if no other d3 is loaded
// This allows microscope pages to load D3 v3 from CDN without conflict
window.d3v7 = d3
if (!window.d3) {
  window.d3 = d3
}

// Export noUiSlider for ichart.js which uses it
export { FPPropChart, LineageChart, noUiSlider }
export default { FPPropChart, LineageChart, noUiSlider }

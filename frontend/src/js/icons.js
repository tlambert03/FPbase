/**
 * JavaScript icon helper for generating SVG icons dynamically
 * This complements the Django {% icon %} template tag for JavaScript-generated HTML
 *
 * Icons are imported from backend/fpbase/static/icons/ using Vite's ?raw import.
 * This ensures single source of truth - same SVGs used in templates and JS.
 * Icons are embedded at build time, no HTTP requests needed.
 */

import alertSvg from "../../../backend/fpbase/static/icons/alert.svg?raw"
import removeSvg from "../../../backend/fpbase/static/icons/remove.svg?raw"
import viewSvg from "../../../backend/fpbase/static/icons/view.svg?raw"

const ICONS = {
  alert: alertSvg,
  view: viewSvg,
  remove: removeSvg,
}

/**
 * Generate an icon SVG element with classes
 * @param {string} name - Icon name (e.g., 'alert', 'view', 'remove')
 * @param {string} classes - CSS classes to apply (e.g., 'me-2 text-danger')
 * @returns {string} SVG HTML string
 */
export function icon(name, classes = "") {
  const svg = ICONS[name]
  if (!svg) {
    console.error(`Icon '${name}' not found. Available icons: ${Object.keys(ICONS).join(", ")}`)
    return ""
  }

  // Add default svg-inline-icon class plus any additional classes
  const allClasses = `svg-inline-icon ${classes}`.trim()

  // Insert classes into the SVG tag
  return svg.replace("<svg", `<svg class="${allClasses}"`)
}

/**
 * Escape HTML to prevent XSS when used in string interpolation
 * @param {string} str - String to escape
 * @returns {string} Escaped string
 */
export function escapeHtml(str) {
  const div = document.createElement("div")
  div.textContent = str
  return div.innerHTML
}

/**
 * Global autocomplete initialization
 *
 * This module provides utilities for automatically initializing
 * autocomplete selects across the application.
 *
 * Usage:
 *   import { initAutocompleteSelects } from './lib/autocomplete-init.js'
 *   initAutocompleteSelects() // Initialize all on page
 *
 * Or for dynamic content:
 *   initAutocompleteSelects(containerElement) // Initialize within container
 */

import { createAutocompleteSelectFromDataset } from "./autocomplete-select.js"

/**
 * Store for tracking initialized selects to avoid double-initialization
 */
const initializedSelects = new WeakSet()

/**
 * Initialize all autocomplete selects within a container
 *
 * Looks for any <select> elements with a data-search-url attribute
 * and enhances them with autocomplete functionality.
 *
 * @param {Element} [container=document] - Container to search within
 * @returns {number} Number of selects initialized
 */
export function initAutocompleteSelects(container = document) {
  const selects = container.querySelectorAll("select[data-search-url]")
  let count = 0

  selects.forEach((select) => {
    // Skip if already initialized
    if (initializedSelects.has(select)) {
      return
    }

    try {
      createAutocompleteSelectFromDataset(select)
      initializedSelects.add(select)
      count++
    } catch (error) {
      console.error("Failed to initialize autocomplete select:", select, error)
    }
  })

  return count
}

/**
 * Initialize autocomplete selects when DOM is ready
 *
 * Call this in your main entry point to automatically initialize
 * all autocomplete selects on page load.
 */
export function initOnReady() {
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", () => {
      initAutocompleteSelects()
    })
  } else {
    initAutocompleteSelects()
  }
}

/**
 * Create a MutationObserver to auto-initialize dynamically added selects
 *
 * Useful for single-page apps or forms that add fields dynamically.
 *
 * @param {Element} [container=document.body] - Container to observe
 * @returns {MutationObserver} The observer instance (call disconnect() to stop)
 */
export function observeForAutocompleteSelects(container = document.body) {
  const observer = new MutationObserver((mutations) => {
    for (const mutation of mutations) {
      for (const node of mutation.addedNodes) {
        if (node.nodeType === Node.ELEMENT_NODE) {
          // Check if the node itself is an autocomplete select
          if (node.matches?.("select[data-search-url]")) {
            initAutocompleteSelects(node.parentElement)
          }
          // Check for autocomplete selects within the node
          else if (node.querySelector) {
            initAutocompleteSelects(node)
          }
        }
      }
    }
  })

  observer.observe(container, {
    childList: true,
    subtree: true,
  })

  return observer
}

import TomSelect from "tom-select"

/**
 * Default configuration for autocomplete selects
 */
const DEFAULT_CONFIG = {
  valueField: "id",
  labelField: "text",
  searchFields: ["text"],
  placeholder: "Type to search...",
  minQueryLength: 2,
  loadThrottle: 300,
  maxOptions: 20,
  noResultsText: "No results found",
  clearButtonTitle: "Clear selection",
}

/**
 * Creates an autocomplete select using TomSelect
 *
 * This provides a consistent interface for creating searchable selects
 * across the application. The implementation uses TomSelect but could
 * be swapped for another library by changing this module.
 *
 * @param {HTMLSelectElement} element - The select element to enhance
 * @param {Object} config - Configuration options
 * @param {string} config.searchUrl - API endpoint for searching (required)
 * @param {string} [config.valueField='id'] - Field name for option values
 * @param {string} [config.labelField='text'] - Field name for option labels
 * @param {string[]} [config.searchFields=['text']] - Fields to search on
 * @param {string} [config.placeholder] - Placeholder text
 * @param {number} [config.minQueryLength=2] - Minimum query length before search
 * @param {number} [config.loadThrottle=300] - Debounce time for searches (ms)
 * @param {number} [config.maxOptions=20] - Maximum number of options to show
 * @param {string} [config.noResultsText] - Text to show when no results found
 * @param {string} [config.clearButtonTitle] - Title for clear button
 * @param {Function} [config.onLoad] - Optional callback when options load
 * @param {Function} [config.formatOption] - Optional custom option formatter
 * @param {Function} [config.formatItem] - Optional custom selected item formatter
 * @returns {TomSelect} The TomSelect instance
 */
export function createAutocompleteSelect(element, config) {
  const {
    searchUrl,
    valueField,
    labelField,
    searchFields,
    placeholder,
    minQueryLength,
    loadThrottle,
    maxOptions,
    noResultsText,
    clearButtonTitle,
    onLoad,
    formatOption,
    formatItem,
  } = { ...DEFAULT_CONFIG, ...config }

  if (!searchUrl) {
    throw new Error("searchUrl is required for autocomplete select")
  }

  return new TomSelect(element, {
    valueField,
    labelField,
    searchField: searchFields,
    placeholder,
    loadThrottle,
    maxOptions,
    plugins: {
      clear_button: {
        title: clearButtonTitle,
      },
    },
    load: (query, callback) => {
      if (!query || query.length < minQueryLength) {
        return callback()
      }

      fetch(`${searchUrl}?q=${encodeURIComponent(query)}`, {
        headers: { "X-Requested-With": "XMLHttpRequest" },
      })
        .then((response) => {
          if (!response.ok) {
            throw new Error(`HTTP error ${response.status}`)
          }
          return response.json()
        })
        .then((data) => {
          const results = data.results || []
          if (onLoad) {
            onLoad(results)
          }
          callback(results)
        })
        .catch(() => {
          callback()
        })
    },
    render: {
      option: formatOption || ((item, escape) => {
        return `<div class="option">${escape(item[labelField])}</div>`
      }),
      item: formatItem || ((item, escape) => {
        return `<div class="item">${escape(item[labelField])}</div>`
      }),
      no_results: () => `<div class="no-results">${noResultsText}</div>`,
    },
  })
}

/**
 * Creates an autocomplete select from data attributes
 *
 * Useful for Alpine.js integration where configuration comes from
 * the DOM rather than JS. Reads config from data-* attributes:
 *
 * - data-search-url (required)
 * - data-value-field
 * - data-label-field
 * - data-search-fields (JSON array)
 * - data-placeholder
 * - data-min-query-length
 *
 * @param {HTMLSelectElement} element - The select element
 * @returns {TomSelect} The TomSelect instance
 */
export function createAutocompleteSelectFromDataset(element) {
  const dataset = element.dataset

  const config = {
    searchUrl: dataset.searchUrl,
  }

  // Optional overrides from data attributes
  if (dataset.valueField) config.valueField = dataset.valueField
  if (dataset.labelField) config.labelField = dataset.labelField
  if (dataset.placeholder) config.placeholder = dataset.placeholder
  if (dataset.minQueryLength) config.minQueryLength = Number.parseInt(dataset.minQueryLength)

  // Parse JSON array for search fields
  if (dataset.searchFields) {
    try {
      config.searchFields = JSON.parse(dataset.searchFields)
    } catch (e) {
      console.warn("Invalid JSON in data-search-fields:", e)
    }
  }

  return createAutocompleteSelect(element, config)
}

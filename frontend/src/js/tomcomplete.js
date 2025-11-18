// src/autocomplete.js
import TomSelect from "tom-select"
import "tom-select/dist/css/tom-select.bootstrap4.css"

export function initAutocomplete(selectElement) {
  const url = selectElement.dataset.autocompleteUrl

  const ts = new TomSelect(selectElement, {
    valueField: "id",
    labelField: "text",
    searchField: "text",

    load: (query, callback) => {
      if (!query.length) return callback()

      fetch(`${url}?q=${encodeURIComponent(query)}`)
        .then((response) => response.json())
        .then((json) => callback(json.results))
        .catch(() => callback())
    },
  })

  // 1Password-specific attribute
  ts.control_input.setAttribute("data-1p-ignore", "")
}

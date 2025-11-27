/**
 * Column Picker UI Module
 *
 * Renders a table showing parsed CSV data and allows users to select
 * which columns contain wavelength and spectrum data.
 */

/**
 * Render the column picker UI.
 *
 * @param {HTMLElement} container - Container element
 * @param {{ headers: string[], rows: number[][] }} parsedCSV - Parsed CSV data
 * @param {function} onColumnsSelected - Callback when columns are selected
 *        Signature: (wavelengthColIndex, dataColIndices[]) => void
 */
export function renderColumnPicker(container, parsedCSV, onColumnsSelected) {
  const state = {
    wavelengthCol: null,
    dataCols: [], // Array of indices for multi-select
  }

  // Clear container
  container.innerHTML = ""
  container.style.display = "block"

  // Instructions
  const instructions = document.createElement("div")
  instructions.className = "alert alert-info mb-3"
  instructions.innerHTML = `
    <strong>Step 2: Select columns</strong><br>
    <span class="text-primary"><i class="bi bi-1-circle me-1"></i>Click a column header to mark it as <strong>Wavelength</strong></span><br>
    <span class="text-info"><i class="bi bi-2-circle me-1"></i>Then click one or more columns to mark them as <strong>Spectrum Data</strong></span>
  `
  container.appendChild(instructions)

  // Status badges
  const statusDiv = document.createElement("div")
  statusDiv.className = "mb-3"
  statusDiv.innerHTML = `
    <span class="badge bg-info me-2" id="wave-status">Wavelength: Not selected</span>
    <span id="data-badges"><span class="badge bg-info">Data: Not selected</span></span>
  `
  container.appendChild(statusDiv)

  // Table wrapper
  const tableWrapper = document.createElement("div")
  tableWrapper.className = "table-responsive"
  tableWrapper.style.maxHeight = "400px"
  tableWrapper.style.overflow = "auto"

  // Create table
  const table = document.createElement("table")
  table.className = "table table-sm table-bordered table-hover"

  // Header row
  const thead = document.createElement("thead")
  thead.className = "sticky-top bg-light"
  const headerRow = document.createElement("tr")

  parsedCSV.headers.forEach((header, index) => {
    const th = document.createElement("th")
    th.scope = "col"
    th.className = "column-header text-center"
    th.style.cursor = "pointer"
    th.style.minWidth = "100px"
    th.innerHTML = `
      <div class="fw-bold">${escapeHtml(header)}</div>
      <small class="text-muted">Col ${index + 1}</small>
    `

    th.addEventListener("click", () => {
      handleColumnClick(index, state, parsedCSV.headers, table, onColumnsSelected)
      updateStatusBadges(state, parsedCSV.headers)
    })

    headerRow.appendChild(th)
  })

  thead.appendChild(headerRow)
  table.appendChild(thead)

  // Data rows (show first 10 for preview)
  const tbody = document.createElement("tbody")
  const previewRows = parsedCSV.rows.slice(0, 10)

  previewRows.forEach((row) => {
    const tr = document.createElement("tr")
    row.forEach((cell) => {
      const td = document.createElement("td")
      td.className = "text-end"
      td.textContent = typeof cell === "number" && !Number.isNaN(cell) ? cell.toFixed(4) : "-"
      tr.appendChild(td)
    })
    tbody.appendChild(tr)
  })

  // Add "more rows" indicator
  if (parsedCSV.rows.length > 10) {
    const tr = document.createElement("tr")
    const td = document.createElement("td")
    td.colSpan = parsedCSV.headers.length
    td.className = "text-center text-muted fst-italic"
    td.textContent = `... and ${parsedCSV.rows.length - 10} more rows`
    tr.appendChild(td)
    tbody.appendChild(tr)
  }

  table.appendChild(tbody)
  tableWrapper.appendChild(table)
  container.appendChild(tableWrapper)

  // Continue button (hidden until valid selection)
  const buttonDiv = document.createElement("div")
  buttonDiv.className = "mt-3"
  buttonDiv.innerHTML = `
    <button type="button" class="btn btn-primary" id="continue-btn" disabled>
      Continue with Selected Columns
    </button>
  `
  container.appendChild(buttonDiv)

  const continueBtn = document.getElementById("continue-btn")
  continueBtn.addEventListener("click", () => {
    if (state.wavelengthCol !== null && state.dataCols.length > 0) {
      onColumnsSelected(state.wavelengthCol, state.dataCols)
    }
  })
}

/**
 * Handle column header click.
 */
function handleColumnClick(index, state, headers, table, _onColumnsSelected) {
  const allThs = table.querySelectorAll("thead th")
  const allTds = table.querySelectorAll("tbody td")
  const columnCount = headers.length

  if (state.wavelengthCol === null) {
    // First click: select wavelength column
    state.wavelengthCol = index
    highlightColumn(allThs, allTds, index, columnCount, "wavelength")
  } else if (state.wavelengthCol === index) {
    // Clicking on wavelength column again: deselect it
    state.wavelengthCol = null
    state.dataCols = []
    clearAllHighlights(allThs, allTds)
  } else if (state.dataCols.includes(index)) {
    // Clicking on already selected data column: deselect it
    state.dataCols = state.dataCols.filter((i) => i !== index)
    clearColumnHighlight(allThs, allTds, index, columnCount)
  } else {
    // Click on a new column: add as data column
    state.dataCols.push(index)
    highlightColumn(allThs, allTds, index, columnCount, "data")
  }

  // Update continue button state
  const continueBtn = document.getElementById("continue-btn")
  continueBtn.disabled = !(state.wavelengthCol !== null && state.dataCols.length > 0)
}

/**
 * Highlight a column with the specified type.
 */
function highlightColumn(allThs, allTds, colIndex, columnCount, type) {
  const className = type === "wavelength" ? "table-primary" : "table-info"

  // Highlight header
  allThs[colIndex].classList.add(className)

  // Highlight all cells in column
  allTds.forEach((td, i) => {
    if (i % columnCount === colIndex) {
      td.classList.add(className)
    }
  })
}

/**
 * Clear highlight from a specific column.
 */
function clearColumnHighlight(allThs, allTds, colIndex, columnCount) {
  allThs[colIndex].classList.remove("table-primary", "table-info")
  allTds.forEach((td, i) => {
    if (i % columnCount === colIndex) {
      td.classList.remove("table-primary", "table-info")
    }
  })
}

/**
 * Clear all column highlights.
 */
function clearAllHighlights(allThs, allTds) {
  allThs.forEach((th) => {
    th.classList.remove("table-primary", "table-info")
  })
  allTds.forEach((td) => {
    td.classList.remove("table-primary", "table-info")
  })
}

/**
 * Update the status badges.
 */
function updateStatusBadges(state, headers) {
  const waveStatus = document.getElementById("wave-status")
  const _dataStatus = document.getElementById("data-status")

  if (state.wavelengthCol !== null) {
    waveStatus.textContent = `Wavelength: ${headers[state.wavelengthCol]}`
    waveStatus.classList.remove("bg-info")
    waveStatus.classList.add("bg-primary")
  } else {
    waveStatus.textContent = "Wavelength: Not selected"
    waveStatus.classList.remove("bg-primary")
    waveStatus.classList.add("bg-info")
  }

  const dataBadges = document.getElementById("data-badges")
  if (state.dataCols.length > 0) {
    // Create individual badge for each data column
    const badges = state.dataCols
      .map((i) => `<span class="badge bg-secondary me-1">${escapeHtml(headers[i])}</span>`)
      .join("")
    dataBadges.innerHTML = badges
  } else {
    dataBadges.innerHTML = `<span class="badge bg-info">Data: Not selected</span>`
  }
}

/**
 * Escape HTML to prevent XSS.
 */
function escapeHtml(str) {
  const div = document.createElement("div")
  div.textContent = str
  return div.innerHTML
}

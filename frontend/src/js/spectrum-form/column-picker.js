/**
 * Column Picker UI
 *
 * Interactive table for selecting wavelength and data columns from parsed CSV.
 */

const PREVIEW_ROWS = 9

/**
 * Render the column picker UI.
 *
 * @param {HTMLElement} container - Container element
 * @param {{ headers: string[], rows: number[][] }} parsedCSV - Parsed CSV data
 * @param {function} onColumnsSelected - Callback: (wavelengthColIndex, dataColIndices[]) => void
 * @param {Object} [options] - Options
 * @param {number} [options.maxColumns] - Maximum number of data columns that can be selected
 */
export function renderColumnPicker(container, parsedCSV, onColumnsSelected, options = {}) {
  const { maxColumns = Infinity } = options
  const state = { wavelengthCol: null, dataCols: [] }

  container.innerHTML = ""
  container.style.display = "block"

  // Instructions
  container.appendChild(createStep2Instructions())

  // Status badges
  const statusDiv = createStatusBadges()
  container.appendChild(statusDiv)

  // Table
  const { table, updateHighlights } = createPreviewTable(parsedCSV, (colIndex) => {
    handleColumnClick(colIndex, state, maxColumns)
    updateHighlights(state)
    updateStatusBadges(state, parsedCSV.headers, maxColumns)
    continueBtn.disabled = !(state.wavelengthCol !== null && state.dataCols.length > 0)
  })
  container.appendChild(table)

  // Continue button
  const buttonDiv = document.createElement("div")
  buttonDiv.className = "mt-3"
  const continueBtn = document.createElement("button")
  continueBtn.type = "button"
  continueBtn.className = "btn btn-primary"
  continueBtn.id = "continue-btn"
  continueBtn.disabled = true
  continueBtn.textContent = "Continue with Selected Columns"
  continueBtn.addEventListener("click", () => {
    if (state.wavelengthCol !== null && state.dataCols.length > 0) {
      onColumnsSelected(state.wavelengthCol, state.dataCols)
    }
  })
  buttonDiv.appendChild(continueBtn)
  container.appendChild(buttonDiv)
}

function createStep2Instructions() {
  const div = document.createElement("div")
  div.className = "alert alert-info mb-3"
  div.innerHTML = `
    <strong>Step 2: Select columns</strong><br>
    <span class="text-primary"><i class="bi bi-1-circle me-1"></i>Click a column header to mark it as <strong>Wavelength</strong></span><br>
    <span class="text-info"><i class="bi bi-2-circle me-1"></i>Then click one or more columns to mark them as <strong>Spectrum Data</strong></span>
  `
  return div
}

function createStatusBadges() {
  const div = document.createElement("div")
  div.className = "mb-3"
  div.innerHTML = `
    <span class="badge bg-info me-2" id="wave-status">Wavelength: Not selected</span>
    <span id="data-badges"><span class="badge bg-info">Data: Not selected</span></span>
  `
  return div
}

function createPreviewTable(parsedCSV, onColumnClick) {
  const wrapper = document.createElement("div")
  wrapper.className = "table-responsive"
  wrapper.style.maxHeight = "420px"
  wrapper.style.overflow = "auto"

  const table = document.createElement("table")
  table.className = "table table-sm table-bordered table-hover"

  // Header
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
    th.addEventListener("click", () => onColumnClick(index))
    headerRow.appendChild(th)
  })

  thead.appendChild(headerRow)
  table.appendChild(thead)

  // Body
  const tbody = document.createElement("tbody")
  const previewRows = parsedCSV.rows.slice(0, PREVIEW_ROWS)

  for (const row of previewRows) {
    const tr = document.createElement("tr")
    for (const cell of row) {
      const td = document.createElement("td")
      td.className = "text-end"
      td.textContent = typeof cell === "number" && !Number.isNaN(cell) ? cell.toFixed(4) : "-"
      tr.appendChild(td)
    }
    tbody.appendChild(tr)
  }

  if (parsedCSV.rows.length > PREVIEW_ROWS) {
    const tr = document.createElement("tr")
    const td = document.createElement("td")
    td.colSpan = parsedCSV.headers.length
    td.className = "text-muted fst-italic"
    td.textContent = `... and ${parsedCSV.rows.length - PREVIEW_ROWS} more rows`
    tr.appendChild(td)
    tbody.appendChild(tr)
  }

  table.appendChild(tbody)
  wrapper.appendChild(table)

  const updateHighlights = (state) => {
    const allThs = table.querySelectorAll("thead th")
    const allTds = table.querySelectorAll("tbody td")
    const columnCount = parsedCSV.headers.length

    // Clear all highlights
    for (const th of allThs) {
      th.classList.remove("table-primary", "table-info")
    }
    for (const td of allTds) {
      td.classList.remove("table-primary", "table-info")
    }

    // Apply wavelength highlight
    if (state.wavelengthCol !== null) {
      highlightColumn(allThs, allTds, state.wavelengthCol, columnCount, "table-primary")
    }

    // Apply data column highlights
    for (const colIndex of state.dataCols) {
      highlightColumn(allThs, allTds, colIndex, columnCount, "table-info")
    }
  }

  return { table: wrapper, updateHighlights }
}

function handleColumnClick(index, state, maxColumns) {
  if (state.wavelengthCol === null) {
    // First click: set wavelength column
    state.wavelengthCol = index
  } else if (state.wavelengthCol === index) {
    // Click on wavelength: deselect everything
    state.wavelengthCol = null
    state.dataCols = []
  } else if (state.dataCols.includes(index)) {
    // Click on selected data column: deselect it
    state.dataCols = state.dataCols.filter((i) => i !== index)
  } else if (state.dataCols.length < maxColumns) {
    // Click on new column: add as data column (if under limit)
    state.dataCols.push(index)
  }
  // If at maxColumns limit, clicking a new column does nothing
}

function highlightColumn(allThs, allTds, colIndex, columnCount, className) {
  allThs[colIndex].classList.add(className)
  allTds.forEach((td, i) => {
    if (i % columnCount === colIndex) {
      td.classList.add(className)
    }
  })
}

function updateStatusBadges(state, headers, maxColumns) {
  const waveStatus = document.getElementById("wave-status")
  const dataBadges = document.getElementById("data-badges")

  if (state.wavelengthCol !== null) {
    waveStatus.textContent = `Wavelength: ${headers[state.wavelengthCol]}`
    waveStatus.className = "badge bg-primary me-2"
  } else {
    waveStatus.textContent = "Wavelength: Not selected"
    waveStatus.className = "badge bg-info me-2"
  }

  if (state.dataCols.length > 0) {
    const limitText =
      maxColumns < Infinity
        ? ` <small class="opacity-75">(${state.dataCols.length}/${maxColumns})</small>`
        : ""
    dataBadges.innerHTML =
      state.dataCols
        .map((i) => `<span class="badge bg-secondary me-1">${escapeHtml(headers[i])}</span>`)
        .join("") + limitText
  } else {
    dataBadges.innerHTML = `<span class="badge bg-info">Data: Not selected</span>`
  }
}

function escapeHtml(str) {
  const div = document.createElement("div")
  div.textContent = str
  return div.innerHTML
}

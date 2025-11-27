/**
 * Spectrum Form Controller
 *
 * Orchestrates the multi-spectrum submission workflow:
 * file upload → column selection → per-spectrum metadata → submission
 */

import { renderColumnPicker } from "./column-picker.js"
import { extractSpectrum, parseCSV } from "./csv-parser.js"
import { interpolateToOneNm, normalize2P, normalizeSpectrum } from "./normalization.js"
import { createSpectrumChart } from "./spectrum-chart.js"

// ============================================================================
// Configuration
// ============================================================================

const CONFIG = {
  peakSearchRadius: 25, // nm around click for peak search
  doiPattern: /^10\.\d{4,}\/[^\s]+$/,
}

/** Categories and their valid subtypes */
const CATEGORY_SUBTYPES = {
  d: [
    { value: "ex", label: "Excitation" },
    { value: "ab", label: "Absorption" },
    { value: "em", label: "Emission" },
    { value: "2p", label: "Two-Photon" },
  ],
  p: [
    { value: "ex", label: "Excitation" },
    { value: "ab", label: "Absorption" },
    { value: "em", label: "Emission" },
    { value: "2p", label: "Two-Photon" },
  ],
  f: [
    { value: "bp", label: "Bandpass" },
    { value: "bx", label: "Bandpass-Ex" },
    { value: "bm", label: "Bandpass-Em" },
    { value: "sp", label: "Shortpass" },
    { value: "lp", label: "Longpass" },
    { value: "bs", label: "Beamsplitter" },
  ],
  c: [{ value: "qe", label: "Quantum Efficiency" }],
  l: [{ value: "pd", label: "Power Distribution" }],
}

const CATEGORY_LABELS = { d: "Dye", p: "Protein", f: "Filter", c: "Camera", l: "Light Source" }

/** Categories that use Select2 autocomplete for owner field */
const AUTOCOMPLETE_URLS = { p: "/autocomplete-protein/" }

/** Categories that show pH/Solvent fields */
const BIO_CATEGORIES = new Set(["d", "p"])

/** Categories that should NOT be normalized (data is percentage-based) */
const NO_NORMALIZE_CATEGORIES = new Set(["f", "c"])

/** Scale factor units by subtype */
const SCALE_UNITS = {
  ex: "EC (M⁻¹cm⁻¹)",
  ab: "EC (M⁻¹cm⁻¹)",
  em: "QE (0-1)",
  "2p": "Cross Section (GM)",
}

// ============================================================================
// Main Entry Point
// ============================================================================

export function initSpectrumForm() {
  const form = document.getElementById("spectrum-form-v2")
  if (!form) return

  const elements = {
    form,
    fileInput: document.getElementById("id_file"),
    columnPicker: document.getElementById("column-picker-container"),
    spectraPreview: document.getElementById("spectra-preview-container"),
    spectraJson: document.getElementById("id_spectra_json"),
    submitBtn: form.querySelector('button[type="submit"]'),
    globalSourceFields: document.getElementById("global-source-fields"),
    confirmationSection: document.getElementById("confirmation-section"),
    sourceInput: document.getElementById("id_source"),
    referenceInput: document.getElementById("id_reference"),
  }

  const state = {
    parsedCSV: null,
    wavelengthCol: null,
    dataCols: [],
    spectra: [],
  }

  setupEventListeners(elements, state)
}

// ============================================================================
// Event Setup
// ============================================================================

function setupEventListeners(el, state) {
  el.referenceInput?.addEventListener("input", (e) => {
    const value = e.target.value.trim()
    const isValid = !value || CONFIG.doiPattern.test(value)
    e.target.setCustomValidity(isValid ? "" : "Please enter a valid DOI (e.g., 10.1234/example)")
    updateFormState(el, state)
  })

  el.sourceInput?.addEventListener("input", () => updateFormState(el, state))

  el.fileInput?.addEventListener("change", async (e) => {
    const file = e.target.files[0]
    if (!file) return

    try {
      const text = await file.text()
      state.parsedCSV = parseCSV(text)

      renderColumnPicker(el.columnPicker, state.parsedCSV, (waveCol, dataCols) => {
        state.wavelengthCol = waveCol
        state.dataCols = dataCols
        processSelectedColumns(el, state)
      })

      // Reset UI until columns are selected
      el.spectraPreview.style.display = "none"
      el.spectraPreview.innerHTML = ""
      el.globalSourceFields.style.display = "none"
      el.confirmationSection.style.display = "none"
      el.submitBtn.disabled = true
    } catch (error) {
      showAlert(el.columnPicker, `Error parsing file: ${error.message}`, "danger")
    }
  })
}

// ============================================================================
// Column Processing
// ============================================================================

function processSelectedColumns(el, state) {
  // Clean up previous spectra
  for (const s of state.spectra) {
    s.chartController?.destroy()
  }
  state.spectra = []

  el.columnPicker.style.display = "none"
  el.spectraPreview.style.display = "block"
  el.spectraPreview.innerHTML = ""

  for (let i = 0; i < state.dataCols.length; i++) {
    const dataColIndex = state.dataCols[i]
    const columnName = state.parsedCSV.headers[dataColIndex]
    const rawData = extractSpectrum(state.parsedCSV, state.wavelengthCol, dataColIndex)

    if (rawData.length < 2) {
      showAlert(el.spectraPreview, `Column "${columnName}" has insufficient data points`, "warning")
      continue
    }

    const spectrum = {
      columnName,
      raw: rawData,
      interpolated: interpolateToOneNm(rawData),
      processed: null,
      chartController: null,
      category: "",
      owner: "",
      subtype: "",
      scaleFactor: null,
      ph: null,
      solvent: "",
      manualPeakWave: null,
    }

    const card = createSpectrumCard(el, state, spectrum, i)
    el.spectraPreview.appendChild(card)
    processSpectrum(spectrum, i)
    state.spectra.push(spectrum)
  }

  el.globalSourceFields.style.display = "block"
  el.confirmationSection.style.display = "block"
  updateFormState(el, state)
}

// ============================================================================
// Spectrum Card Creation
// ============================================================================

function createSpectrumCard(el, state, spectrum, index) {
  const card = document.createElement("div")
  card.className = "card mb-3 spectrum-card"
  card.id = `spectrum-card-${index}`

  const [minWave, maxWave] = [
    Math.round(spectrum.raw[0][0]),
    Math.round(spectrum.raw[spectrum.raw.length - 1][0]),
  ]

  const shouldNormalize = !NO_NORMALIZE_CATEGORIES.has(spectrum.category)
  const showBioFields = BIO_CATEGORIES.has(spectrum.category)
  const useAutocomplete = spectrum.category in AUTOCOMPLETE_URLS

  card.innerHTML = buildCardHTML(spectrum, index, {
    shouldNormalize,
    showBioFields,
    useAutocomplete,
  })

  // Attach event handlers after next tick (DOM must be ready)
  setTimeout(() => attachCardEventHandlers(el, state, spectrum, index), 0)

  // Store card info for footer
  card.dataset.minWave = minWave
  card.dataset.maxWave = maxWave
  card.dataset.pointCount = spectrum.raw.length

  return card
}

function buildCardHTML(spectrum, index, { shouldNormalize, showBioFields, useAutocomplete }) {
  const categoryOptions = buildCategoryOptions(spectrum.category)
  const subtypeOptions = buildSubtypeOptions(spectrum.category, spectrum.subtype)
  const scaleUnits = SCALE_UNITS[spectrum.subtype] || ""
  const [minWave, maxWave] = [
    Math.round(spectrum.raw[0][0]),
    Math.round(spectrum.raw[spectrum.raw.length - 1][0]),
  ]

  return `
    <div class="card-header d-flex justify-content-between align-items-center">
      <strong>Spectrum ${index + 1}: ${escapeHtml(spectrum.columnName)}</strong>
      <span class="badge bg-secondary" id="peak-badge-${index}"
            style="${shouldNormalize ? "" : "display: none;"}">Peak: --</span>
    </div>
    <div class="card-body">
      <div class="row mb-3">
        <div class="col-md-6">
          <label class="form-label" id="category-label-${index}">
            Spectrum Category <span class="text-danger">*</span>
          </label>
          <select class="form-select form-select-sm" id="category-select-${index}" required>
            ${categoryOptions}
          </select>
        </div>
        <div class="col-md-6">
          <label class="form-label" id="subtype-label-${index}">
            Subtype <span class="text-danger">*</span>
          </label>
          <select class="form-select form-select-sm" id="subtype-select-${index}" required>
            ${subtypeOptions}
          </select>
        </div>
      </div>

      <div class="row mb-3">
        <div class="col-12">
          <label class="form-label" id="owner-label-${index}">
            Owner <span class="text-danger">*</span>
          </label>
          <input type="text" class="form-control form-control-sm" id="owner-input-${index}"
                 placeholder="Name of dye, filter, etc." required
                 style="${useAutocomplete ? "display: none;" : ""}">
          <select class="form-select form-select-sm" id="owner-select-${index}"
                  style="${useAutocomplete ? "" : "display: none;"}">
            <option value="">Select a protein...</option>
          </select>
        </div>
      </div>

      <div class="row mb-3 optional-fields-${index}"
           style="${shouldNormalize || showBioFields ? "" : "display: none;"}">
        <div class="col-md-4" id="scale-factor-container-${index}"
             style="${shouldNormalize ? "" : "display: none;"}">
          <label class="form-label">Scale Factor <small class="text-muted">(optional)</small></label>
          <div class="d-flex align-items-center gap-2">
            <input type="number" class="form-control form-control-sm no-spinners flex-grow-1"
                   id="scale-factor-${index}" step="any">
            <small class="text-muted text-nowrap flex-shrink-0"
                   id="scale-factor-units-${index}">${scaleUnits}</small>
          </div>
          <div class="form-text small">Absolute magnitude at peak wavelength</div>
        </div>
        <div class="col-md-4 bio-field-${index}" style="${showBioFields ? "" : "display: none;"}">
          <label class="form-label">pH <small class="text-muted">(optional)</small></label>
          <input type="number" class="form-control form-control-sm" id="ph-input-${index}"
                 step="0.1" min="0" max="14" placeholder="e.g., 7.4">
        </div>
        <div class="col-md-4 bio-field-${index}" style="${showBioFields ? "" : "display: none;"}">
          <label class="form-label">Solvent <small class="text-muted">(optional)</small></label>
          <input type="text" class="form-control form-control-sm" id="solvent-input-${index}"
                 placeholder="e.g., PBS, DMSO">
        </div>
      </div>

      <div id="chart-container-${index}"
           style="cursor: ${shouldNormalize ? "crosshair" : "default"};"></div>
      <div class="small text-muted" id="chart-hint-${index}"
           style="${shouldNormalize ? "" : "display: none;"}">
        <em>Click on the chart to set peak location (searches ±${CONFIG.peakSearchRadius}nm)</em>
      </div>

      <div class="mt-2 d-flex justify-content-between align-items-center">
        <span class="text-muted small">${spectrum.raw.length} points | ${minWave}-${maxWave} nm</span>
        <button type="button" class="btn btn-sm btn-outline-danger" id="remove-btn-${index}">
          Remove
        </button>
      </div>
    </div>
  `
}

function attachCardEventHandlers(el, state, spectrum, index) {
  const categorySelect = document.getElementById(`category-select-${index}`)
  const subtypeSelect = document.getElementById(`subtype-select-${index}`)
  const ownerInput = document.getElementById(`owner-input-${index}`)
  const ownerSelect = document.getElementById(`owner-select-${index}`)
  const scaleFactorInput = document.getElementById(`scale-factor-${index}`)
  const phInput = document.getElementById(`ph-input-${index}`)
  const solventInput = document.getElementById(`solvent-input-${index}`)
  const removeBtn = document.getElementById(`remove-btn-${index}`)

  // Select2 initialization for protein autocomplete
  const initOwnerSelect2 = (category) => {
    const url = AUTOCOMPLETE_URLS[category]
    if (!url || !window.$ || !ownerSelect) return

    if ($(ownerSelect).hasClass("select2-hidden-accessible")) {
      $(ownerSelect).select2("destroy")
    }

    $(ownerSelect).select2({
      theme: "bootstrap-5",
      placeholder: "Search for a protein...",
      allowClear: true,
      ajax: {
        url,
        dataType: "json",
        delay: 250,
        cache: true,
        data: (params) => ({ q: params.term, page: params.page }),
        processResults: (data) => data,
      },
    })

    $(ownerSelect).on("select2:select select2:clear", () => {
      const selected = $(ownerSelect).select2("data")
      spectrum.owner = selected?.[0]?.text || ""
      updateFormState(el, state)
    })
  }

  if (spectrum.category in AUTOCOMPLETE_URLS) {
    initOwnerSelect2(spectrum.category)
  }

  categorySelect?.addEventListener("change", (e) => {
    const newCategory = e.target.value
    spectrum.category = newCategory

    // Update subtypes (auto-select if only one option)
    const subtypes = CATEGORY_SUBTYPES[newCategory] || []
    spectrum.subtype = subtypes.length === 1 ? subtypes[0].value : ""
    subtypeSelect.innerHTML = buildSubtypeOptions(newCategory, spectrum.subtype)

    const shouldNorm = !NO_NORMALIZE_CATEGORIES.has(newCategory)
    const showBio = BIO_CATEGORIES.has(newCategory)
    const useAutocomplete = newCategory in AUTOCOMPLETE_URLS

    // Toggle owner input/select visibility
    if (ownerInput) ownerInput.style.display = useAutocomplete ? "none" : ""
    if (ownerSelect) ownerSelect.style.display = useAutocomplete ? "" : "none"

    // Reset owner when switching modes
    spectrum.owner = ""
    if (ownerInput) ownerInput.value = ""
    if (ownerSelect && $(ownerSelect).hasClass("select2-hidden-accessible")) {
      $(ownerSelect).val(null).trigger("change")
    }

    if (useAutocomplete) {
      initOwnerSelect2(newCategory)
    } else if (ownerSelect && window.$ && $(ownerSelect).hasClass("select2-hidden-accessible")) {
      $(ownerSelect).select2("destroy")
    }

    // Update UI visibility
    updateCardVisibility(index, { shouldNorm, showBio })
    updateScaleFactorUnits(index, spectrum.subtype)

    spectrum.manualPeakWave = null
    processSpectrum(spectrum, index)
    updateFormState(el, state)
  })

  subtypeSelect?.addEventListener("change", (e) => {
    spectrum.subtype = e.target.value
    spectrum.manualPeakWave = null
    updateScaleFactorUnits(index, spectrum.subtype)
    processSpectrum(spectrum, index)
    updateFormState(el, state)
  })

  ownerInput?.addEventListener("input", (e) => {
    spectrum.owner = e.target.value.trim()
    updateFormState(el, state)
  })

  scaleFactorInput?.addEventListener("change", (e) => {
    const val = parseFloat(e.target.value)
    spectrum.scaleFactor = Number.isNaN(val) ? null : val
    updateFormState(el, state)
  })

  phInput?.addEventListener("change", (e) => {
    const val = parseFloat(e.target.value)
    spectrum.ph = Number.isNaN(val) ? null : val
    updateFormState(el, state)
  })

  solventInput?.addEventListener("input", (e) => {
    spectrum.solvent = e.target.value.trim()
    updateFormState(el, state)
  })

  removeBtn?.addEventListener("click", () => {
    if (ownerSelect && window.$ && $(ownerSelect).hasClass("select2-hidden-accessible")) {
      $(ownerSelect).select2("destroy")
    }
    spectrum.chartController?.destroy()
    document.getElementById(`spectrum-card-${index}`)?.remove()
    state.spectra = state.spectra.filter((s) => s !== spectrum)
    updateFormState(el, state)

    if (state.spectra.length === 0) {
      el.globalSourceFields.style.display = "none"
      el.confirmationSection.style.display = "none"
    }
  })
}

function updateCardVisibility(index, { shouldNorm, showBio }) {
  const peakBadge = document.getElementById(`peak-badge-${index}`)
  const scaleContainer = document.getElementById(`scale-factor-container-${index}`)
  const optionalRow = document.querySelector(`.optional-fields-${index}`)
  const chartContainer = document.getElementById(`chart-container-${index}`)
  const chartHint = document.getElementById(`chart-hint-${index}`)

  if (peakBadge) peakBadge.style.display = shouldNorm ? "" : "none"
  if (scaleContainer) scaleContainer.style.display = shouldNorm ? "" : "none"
  if (optionalRow) optionalRow.style.display = shouldNorm || showBio ? "" : "none"
  if (chartContainer) chartContainer.style.cursor = shouldNorm ? "crosshair" : "default"
  if (chartHint) chartHint.style.display = shouldNorm ? "" : "none"

  document.querySelectorAll(`.bio-field-${index}`).forEach((field) => {
    field.style.display = showBio ? "" : "none"
  })
}

function updateScaleFactorUnits(index, subtype) {
  const unitsEl = document.getElementById(`scale-factor-units-${index}`)
  if (unitsEl) unitsEl.textContent = SCALE_UNITS[subtype] || ""
}

// ============================================================================
// Spectrum Processing
// ============================================================================

function processSpectrum(spectrum, index) {
  const shouldNormalize = !NO_NORMALIZE_CATEGORIES.has(spectrum.category)
  let processedData
  let peakWave = null

  if (shouldNormalize) {
    const options =
      spectrum.manualPeakWave != null
        ? {
            rangeMin: spectrum.manualPeakWave - CONFIG.peakSearchRadius,
            rangeMax: spectrum.manualPeakWave + CONFIG.peakSearchRadius,
          }
        : {}

    const normFn = spectrum.subtype === "2p" ? normalize2P : normalizeSpectrum
    const result = normFn(spectrum.interpolated, options)
    processedData = result.normalized
    peakWave = result.peakWave
  } else {
    // For filters/cameras: convert percentage if needed
    const maxVal = Math.max(...spectrum.interpolated.map(([, y]) => y))
    processedData =
      maxVal > 1 ? spectrum.interpolated.map(([x, y]) => [x, y / 100]) : spectrum.interpolated
  }

  spectrum.processed = processedData

  // Update peak badge
  const peakBadge = document.getElementById(`peak-badge-${index}`)
  if (peakBadge) {
    peakBadge.textContent = peakWave !== null ? `Peak: ${peakWave} nm` : "No peak"
    peakBadge.className = `badge ${peakWave !== null ? "bg-primary" : "bg-secondary"}`
  }

  // Create or update chart
  const chartContainer = document.getElementById(`chart-container-${index}`)
  const chartOptions = {
    name: spectrum.columnName,
    onClick: shouldNormalize
      ? (wavelength) => {
          spectrum.manualPeakWave = Math.round(wavelength)
          processSpectrum(spectrum, index)
        }
      : null,
  }

  if (!spectrum.chartController && chartContainer) {
    spectrum.chartController = createSpectrumChart(chartContainer, processedData, chartOptions)
  } else if (spectrum.chartController) {
    spectrum.chartController.updateData(processedData, spectrum.columnName)
  }

  if (spectrum.chartController) {
    if (peakWave !== null) {
      spectrum.chartController.setPeakMarker(peakWave)
    } else {
      spectrum.chartController.clearAnnotations()
    }
  }
}

// ============================================================================
// Form State Management
// ============================================================================

function updateFormState(el, state) {
  const validSpectra = state.spectra.filter(
    (s) => s.processed?.length > 0 && s.category && s.subtype && s.owner.trim()
  )

  const spectraJson = validSpectra.map((s) => ({
    data: s.processed,
    category: s.category,
    owner: s.owner,
    subtype: s.subtype,
    scale_factor: s.scaleFactor,
    ph: BIO_CATEGORIES.has(s.category) ? s.ph : null,
    solvent: BIO_CATEGORIES.has(s.category) ? s.solvent : null,
    peak_wave: getPeakWave(s.processed),
    column_name: s.columnName,
  }))

  el.spectraJson.value = JSON.stringify(spectraJson)

  // Track missing fields per spectrum
  const missing = { owners: [], categories: [], subtypes: [] }
  for (const [i, s] of state.spectra.entries()) {
    const hasOwner = s.owner.trim() !== ""
    const hasCategory = s.category !== ""
    const hasSubtype = s.subtype !== ""
    const isComplete = hasOwner && hasCategory && hasSubtype

    updateStatusIcon(i, isComplete)
    updateFieldLabels(i, { hasOwner, hasCategory, hasSubtype })

    if (!hasCategory) missing.categories.push(s.columnName)
    else if (!hasSubtype) missing.subtypes.push(s.columnName)
    if (!hasOwner) missing.owners.push(s.columnName)
  }

  // Check source/reference validation
  const hasSource = el.sourceInput?.value.trim() !== ""
  const hasReference = el.referenceInput?.value.trim() !== ""
  const hasValidReference = !hasReference || CONFIG.doiPattern.test(el.referenceInput.value.trim())
  const hasSourceOrRef = hasSource || hasReference

  // Update submit button
  const allComplete =
    missing.owners.length === 0 && missing.categories.length === 0 && missing.subtypes.length === 0
  const isValid = state.spectra.length > 0 && allComplete && hasSourceOrRef && hasValidReference
  el.submitBtn.disabled = !isValid

  updateValidationMessage(el, state, missing, hasSourceOrRef, hasValidReference)
}

function updateStatusIcon(index, isComplete) {
  const cardHeader = document.querySelector(`#spectrum-card-${index} .card-header`)
  let statusIcon = document.getElementById(`status-icon-${index}`)

  if (!statusIcon && cardHeader) {
    statusIcon = document.createElement("span")
    statusIcon.id = `status-icon-${index}`
    statusIcon.className = "me-2"
    cardHeader.querySelector("strong")?.prepend(statusIcon)
  }

  if (statusIcon) {
    const icon = isComplete ? "✓" : "!"
    const color = isComplete ? "text-success" : "text-warning"
    statusIcon.innerHTML = `<span class="${color}">${icon}</span> `
  }
}

function updateFieldLabels(index, { hasOwner, hasCategory, hasSubtype }) {
  const ownerLabel = document.getElementById(`owner-label-${index}`)
  const categoryLabel = document.getElementById(`category-label-${index}`)
  const subtypeLabel = document.getElementById(`subtype-label-${index}`)

  if (ownerLabel) ownerLabel.style.fontWeight = hasOwner ? "normal" : "bold"
  if (categoryLabel) categoryLabel.style.fontWeight = hasCategory ? "normal" : "bold"
  if (subtypeLabel) subtypeLabel.style.fontWeight = hasSubtype ? "normal" : "bold"
}

function updateValidationMessage(el, state, missing, hasSourceOrRef, hasValidReference) {
  let messageEl = document.getElementById("validation-message")
  if (!messageEl) {
    messageEl = document.createElement("div")
    messageEl.id = "validation-message"
    messageEl.className = "mb-3"
    el.submitBtn.parentElement.parentElement.insertBefore(messageEl, el.submitBtn.parentElement)
  }

  if (state.spectra.length === 0) {
    messageEl.innerHTML = ""
    return
  }

  const issues = []

  if (missing.categories.length > 0) {
    const n = missing.categories.length
    const names = missing.categories.map((name) => `"${escapeHtml(name)}"`).join(", ")
    issues.push(
      `<strong>${n}</strong> ${n === 1 ? "spectrum needs" : "spectra need"} a category: ${names}`
    )
  }

  if (missing.subtypes.length > 0) {
    const n = missing.subtypes.length
    const names = missing.subtypes.map((name) => `"${escapeHtml(name)}"`).join(", ")
    issues.push(
      `<strong>${n}</strong> ${n === 1 ? "spectrum needs" : "spectra need"} a subtype: ${names}`
    )
  }

  if (missing.owners.length > 0) {
    const n = missing.owners.length
    const names = missing.owners.map((name) => `"${escapeHtml(name)}"`).join(", ")
    issues.push(
      `<strong>${n}</strong> ${n === 1 ? "spectrum needs" : "spectra need"} an owner: ${names}`
    )
  }

  if (!hasSourceOrRef) {
    issues.push(
      "Please provide at least one of <strong>Source</strong> or <strong>Primary Reference</strong>"
    )
  } else if (!hasValidReference) {
    issues.push("Primary Reference must be a valid DOI (e.g., 10.1234/example)")
  }

  if (issues.length > 0) {
    messageEl.innerHTML = `
      <div class="alert alert-warning py-2 mb-2">
        <ul class="mb-0 ps-3">${issues.map((i) => `<li>${i}</li>`).join("")}</ul>
      </div>
    `
  } else {
    messageEl.innerHTML = `
      <div class="alert alert-success py-2 mb-2">✓ All spectra ready to submit</div>
    `
  }
}

// ============================================================================
// Utilities
// ============================================================================

function buildCategoryOptions(selectedValue) {
  const placeholder = `<option value=""${!selectedValue ? " selected" : ""}>-----</option>`
  const options = Object.entries(CATEGORY_LABELS)
    .map(
      ([value, label]) =>
        `<option value="${value}"${value === selectedValue ? " selected" : ""}>${label}</option>`
    )
    .join("")
  return placeholder + options
}

function buildSubtypeOptions(category, selectedValue) {
  if (!category) {
    return `<option value="" selected>Select category first...</option>`
  }
  const subtypes = CATEGORY_SUBTYPES[category] || []
  if (subtypes.length === 1) {
    return `<option value="${subtypes[0].value}" selected>${subtypes[0].label}</option>`
  }
  const placeholder = `<option value=""${!selectedValue ? " selected" : ""}>-----</option>`
  const options = subtypes
    .map(
      ({ value, label }) =>
        `<option value="${value}"${value === selectedValue ? " selected" : ""}>${label}</option>`
    )
    .join("")
  return placeholder + options
}

function getPeakWave(data) {
  if (!data?.length) return null
  let maxVal = -Infinity
  let peakWave = null
  for (const [wave, val] of data) {
    if (val > maxVal) {
      maxVal = val
      peakWave = wave
    }
  }
  return peakWave
}

function showAlert(container, message, type = "info") {
  const alert = document.createElement("div")
  alert.className = `alert alert-${type} alert-dismissible fade show`
  alert.innerHTML = `${message}<button type="button" class="btn-close" data-bs-dismiss="alert"></button>`
  container.prepend(alert)
}

function escapeHtml(str) {
  const div = document.createElement("div")
  div.textContent = str
  return div.innerHTML
}

/**
 * Spectrum Form Controller
 *
 * Orchestrates the multi-spectrum submission workflow:
 * file upload → column selection → per-spectrum metadata → submission
 */

import { renderColumnPicker } from "./column-picker.js"
import { extractSpectrum, parseCSV } from "./csv-parser.js"
import { checkSimilarOwners } from "./duplicate-checker.js"
import { interpolateToOneNm, normalize2P, normalizeSpectrum } from "./normalization.js"
import { createSpectrumChart } from "./spectrum-chart.js"

// ============================================================================
// Type Definitions
// ============================================================================

/**
 * @typedef {Object} SpectrumJSON
 * JSON structure sent to the server on form submission.
 * Must match SpectrumJSONData in backend/proteins/forms/spectrum_v2.py.
 * @property {Array<[number, number]>} data - Array of [wavelength, value] pairs
 * @property {string} category - Category code (d=dye, p=protein, f=filter, c=camera, l=light)
 * @property {string} owner - Owner name (protein/dye/filter name)
 * @property {string} subtype - Subtype code (ex, ab, em, 2p, bp, etc.)
 * @property {number|null} scale_factor - Scale factor at peak wavelength
 * @property {number|null} ph - pH value (for fluorophore categories only)
 * @property {string|null} solvent - Solvent name (for fluorophore categories only)
 * @property {number|null} peak_wave - Peak wavelength in nm
 * @property {string} column_name - Original column name from CSV
 */

/**
 * @typedef {Object} Spectrum
 * Internal spectrum object used during form editing.
 * @property {string} columnName - Original column name from CSV
 * @property {Array<[number, number]>} raw - Raw data points from CSV
 * @property {Array<[number, number]>} interpolated - Interpolated to 1nm spacing
 * @property {Array<[number, number]>|null} processed - Normalized/processed data
 * @property {Object|null} chartController - Chart.js controller instance
 * @property {string} category - Category code
 * @property {string} owner - Owner name
 * @property {string} subtype - Subtype code
 * @property {number|null} scaleFactor - Scale factor value
 * @property {number|null} ph - pH value
 * @property {string} solvent - Solvent name
 * @property {number|null} manualPeakWave - User-selected peak wavelength
 * @property {boolean} hasExactMatch - Whether an exact match exists in database
 */

// ============================================================================
// Configuration
// ============================================================================

const CONFIG = {
  peakSearchRadius: 25, // nm around click for peak search
  doiPattern: /^10\.\d{4,}\/[^\s]+$/,
}

/** Subtypes available for fluorophore spectra (proteins and dyes) */
const FLUOR_SUBTYPES = [
  { value: "ex", label: "Excitation" },
  { value: "ab", label: "Absorption" },
  { value: "em", label: "Emission" },
  { value: "2p", label: "Two-Photon" },
]

/** Categories and their valid subtypes */
const CATEGORY_SUBTYPES = {
  d: FLUOR_SUBTYPES,
  p: FLUOR_SUBTYPES,
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

/** Human-readable category labels */
const CATEGORY_LABELS = {
  d: "Dye",
  p: "Protein",
  f: "Filter",
  c: "Camera",
  l: "Light Source",
}

/** Categories that use Select2 autocomplete for owner field */
const AUTOCOMPLETE_URLS = { p: "/autocomplete-protein/" }

/** Categories that show pH/Solvent fields */
const FLUOR_CATEGORIES = new Set(["d", "p"])

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
// Category Helper Functions
// ============================================================================

/**
 * Determine display/behavior flags for a given category.
 * @param {string} category - Category code
 * @returns {{shouldNormalize: boolean, showFluor: boolean, useAutocomplete: boolean}}
 */
function getCategoryFlags(category) {
  return {
    shouldNormalize: !!category && !NO_NORMALIZE_CATEGORIES.has(category),
    showFluor: FLUOR_CATEGORIES.has(category),
    useAutocomplete: category in AUTOCOMPLETE_URLS,
  }
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
    confirmationCheckbox: document.getElementById("id_confirmation"),
  }

  const state = {
    parsedCSV: null,
    wavelengthCol: null,
    dataCols: [],
    spectra: [],
  }

  setupEventListeners(elements, state)

  // Restore state from spectra_json if present (e.g., after form validation error)
  restoreStateFromJson(elements, state)
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
  el.confirmationCheckbox?.addEventListener("change", () => updateFormState(el, state))

  el.fileInput?.addEventListener("change", async (e) => {
    const file = e.target.files[0]
    if (!file) return

    resetFormUI(el, state)

    try {
      const text = await file.text()
      state.parsedCSV = parseCSV(text)

      renderColumnPicker(el.columnPicker, state.parsedCSV, (waveCol, dataCols) => {
        state.wavelengthCol = waveCol
        state.dataCols = dataCols
        processSelectedColumns(el, state)
      })
    } catch (error) {
      el.columnPicker.innerHTML = ""
      el.columnPicker.style.display = "block"
      showAlert(el.columnPicker, `Error reading file: ${error.message}`, "danger")
      el.fileInput.value = ""
    }
  })
}

/**
 * Reset the form UI to initial state (before file selection).
 */
function resetFormUI(el, state) {
  el.spectraPreview.style.display = "none"
  el.spectraPreview.innerHTML = ""
  el.globalSourceFields.style.display = "none"
  el.confirmationSection.style.display = "none"
  el.submitBtn.disabled = true

  document.getElementById("form-errors")?.remove()
  document.getElementById("validation-message")?.remove()

  for (const s of state.spectra) {
    s.chartController?.destroy()
  }
  state.spectra = []
  state.parsedCSV = null
  state.wavelengthCol = null
  state.dataCols = []
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
  el.spectraPreview.appendChild(createStep3Instructions())

  for (let i = 0; i < state.dataCols.length; i++) {
    const dataColIndex = state.dataCols[i]
    const columnName = state.parsedCSV.headers[dataColIndex]
    const rawData = extractSpectrum(state.parsedCSV, state.wavelengthCol, dataColIndex)

    if (rawData.length < 2) {
      showAlert(el.spectraPreview, `Column "${columnName}" has insufficient data points`, "warning")
      continue
    }

    const spectrum = createSpectrumObject(columnName, rawData)
    const card = createSpectrumCard(el, state, spectrum, i)
    el.spectraPreview.appendChild(card)
    processSpectrum(spectrum, i, el, state)
    state.spectra.push(spectrum)
  }

  el.globalSourceFields.style.display = "block"
  el.confirmationSection.style.display = "block"
  updateFormState(el, state)
}

/**
 * Create a new spectrum object with default values.
 * @param {string} columnName - Column name from CSV
 * @param {Array<[number, number]>} rawData - Raw data points
 * @returns {Spectrum}
 */
function createSpectrumObject(columnName, rawData) {
  return {
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
    hasExactMatch: false,
  }
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

  const flags = getCategoryFlags(spectrum.category)
  card.innerHTML = buildCardHTML(spectrum, index, flags)

  // Attach event handlers after next tick (DOM must be ready)
  setTimeout(() => attachCardEventHandlers(el, state, spectrum, index), 0)

  card.dataset.minWave = minWave
  card.dataset.maxWave = maxWave
  card.dataset.pointCount = spectrum.raw.length

  return card
}

function buildCardHTML(spectrum, index, { shouldNormalize, showFluor, useAutocomplete }) {
  const categoryOptions = buildCategoryOptions(spectrum.category)
  const subtypeOptions = buildSubtypeOptions(spectrum.category, spectrum.subtype)
  const scaleUnits = SCALE_UNITS[spectrum.subtype] || ""
  const [minWave, maxWave] = [
    Math.round(spectrum.raw[0][0]),
    Math.round(spectrum.raw[spectrum.raw.length - 1][0]),
  ]

  const hasCategory = !!spectrum.category
  const ownerPlaceholder = hasCategory ? "Name of dye, filter, etc." : "Select a category first"

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
                 placeholder="${ownerPlaceholder}"
                 ${useAutocomplete ? "" : "required"}
                 ${hasCategory ? "" : "disabled"}
                 style="${useAutocomplete ? "display: none;" : ""}">
          <select class="form-select form-select-sm" id="owner-select-${index}"
                  ${useAutocomplete ? "required" : ""}
                  ${hasCategory ? "" : "disabled"}
                  style="${useAutocomplete ? "" : "display: none;"}">
            <option value="">Select a protein...</option>
          </select>
          <div class="alert alert-warning small mt-2 mb-0 d-none" id="owner-warning-${index}"></div>
        </div>
      </div>

      <div class="row mb-3 optional-fields-${index}"
           style="${shouldNormalize || showFluor ? "" : "display: none;"}">
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
        <div class="col-md-4 fluor-field-${index}" style="${showFluor ? "" : "display: none;"}">
          <label class="form-label">pH <small class="text-muted">(optional)</small></label>
          <input type="number" class="form-control form-control-sm" id="ph-input-${index}"
                 step="0.1" min="0" max="14" placeholder="e.g., 7.4">
        </div>
        <div class="col-md-4 fluor-field-${index}" style="${showFluor ? "" : "display: none;"}">
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

// ============================================================================
// Card Event Handlers
// ============================================================================

function attachCardEventHandlers(el, state, spectrum, index) {
  const categorySelect = document.getElementById(`category-select-${index}`)
  const subtypeSelect = document.getElementById(`subtype-select-${index}`)
  const ownerInput = document.getElementById(`owner-input-${index}`)
  const ownerSelect = document.getElementById(`owner-select-${index}`)
  const scaleFactorInput = document.getElementById(`scale-factor-${index}`)
  const phInput = document.getElementById(`ph-input-${index}`)
  const solventInput = document.getElementById(`solvent-input-${index}`)
  const removeBtn = document.getElementById(`remove-btn-${index}`)

  // Restore saved values to form inputs
  restoreCardInputValues(spectrum, index, el, state)

  // Category change handler
  categorySelect?.addEventListener("change", (e) => {
    handleCategoryChange(e.target.value, spectrum, index, el, state, {
      subtypeSelect,
      ownerInput,
      ownerSelect,
    })
  })

  // Subtype change handler
  subtypeSelect?.addEventListener("change", (e) => {
    handleSubtypeChange(e.target.value, spectrum, index, el, state)
  })

  // Owner input handler (text input for non-protein categories)
  ownerInput?.addEventListener("input", (e) => {
    handleOwnerInput(e.target.value, spectrum, index, el, state)
  })

  // Optional field handlers
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

  // Remove button handler
  removeBtn?.addEventListener("click", () => {
    handleRemoveSpectrum(spectrum, index, el, state, ownerSelect)
  })
}

/**
 * Restore saved values to card input fields (used during card creation).
 */
function restoreCardInputValues(spectrum, index, el, state) {
  const ownerInput = document.getElementById(`owner-input-${index}`)
  const ownerSelect = document.getElementById(`owner-select-${index}`)
  const scaleFactorInput = document.getElementById(`scale-factor-${index}`)
  const phInput = document.getElementById(`ph-input-${index}`)
  const solventInput = document.getElementById(`solvent-input-${index}`)

  // Initialize Select2 or populate text input based on category
  if (spectrum.category in AUTOCOMPLETE_URLS) {
    initOwnerSelect2(ownerSelect, spectrum.category, spectrum.owner, spectrum, el, state)
  } else if (ownerInput && spectrum.owner) {
    ownerInput.value = spectrum.owner
  }

  // Populate optional fields
  if (scaleFactorInput && spectrum.scaleFactor != null) {
    scaleFactorInput.value = spectrum.scaleFactor
  }
  if (phInput && spectrum.ph != null) {
    phInput.value = spectrum.ph
  }
  if (solventInput && spectrum.solvent) {
    solventInput.value = spectrum.solvent
  }
}

/**
 * Handle category selection change.
 */
async function handleCategoryChange(newCategory, spectrum, index, el, state, elements) {
  const { subtypeSelect, ownerInput, ownerSelect } = elements

  spectrum.category = newCategory

  // Update subtypes (auto-select if only one option)
  const subtypes = CATEGORY_SUBTYPES[newCategory] || []
  spectrum.subtype = subtypes.length === 1 ? subtypes[0].value : ""
  subtypeSelect.innerHTML = buildSubtypeOptions(newCategory, spectrum.subtype)

  const flags = getCategoryFlags(newCategory)

  // Toggle owner input/select visibility and required attribute
  configureOwnerFields(ownerInput, ownerSelect, flags.useAutocomplete, !!newCategory)

  // Reset owner when switching modes
  spectrum.owner = ""
  if (ownerInput) ownerInput.value = ""
  if (isSelect2Initialized(ownerSelect)) {
    $(ownerSelect).val(null).trigger("change")
  }

  // Initialize or destroy Select2 based on category
  if (flags.useAutocomplete) {
    initOwnerSelect2(ownerSelect, newCategory, null, spectrum, el, state)
  } else {
    destroySelect2(ownerSelect)
  }

  // Update UI visibility
  updateCardVisibility(index, flags)
  updateScaleFactorUnits(index, spectrum.subtype)

  // Clear duplicate warning and re-check
  await updateOwnerWarning(spectrum, index)

  spectrum.manualPeakWave = null
  processSpectrum(spectrum, index, el, state)
  updateFormState(el, state)
}

/**
 * Handle subtype selection change.
 */
async function handleSubtypeChange(newSubtype, spectrum, index, el, state) {
  spectrum.subtype = newSubtype
  spectrum.manualPeakWave = null
  updateScaleFactorUnits(index, newSubtype)

  await updateOwnerWarning(spectrum, index)

  processSpectrum(spectrum, index, el, state)
  updateFormState(el, state)
}

/**
 * Handle owner text input change.
 */
async function handleOwnerInput(value, spectrum, index, el, state) {
  spectrum.owner = value.trim()
  await updateOwnerWarning(spectrum, index)
  updateFormState(el, state)
}

/**
 * Handle spectrum removal.
 */
function handleRemoveSpectrum(spectrum, index, el, state, ownerSelect) {
  destroySelect2(ownerSelect)
  spectrum.chartController?.destroy()
  document.getElementById(`spectrum-card-${index}`)?.remove()
  state.spectra = state.spectra.filter((s) => s !== spectrum)
  updateFormState(el, state)

  if (state.spectra.length === 0) {
    // Reset to column picker view
    el.spectraPreview.style.display = "none"
    el.spectraPreview.innerHTML = ""
    el.globalSourceFields.style.display = "none"
    el.confirmationSection.style.display = "none"
    el.columnPicker.style.display = "block"
    el.fileInput.value = ""
  }
}

// ============================================================================
// Owner Field Helpers
// ============================================================================

/**
 * Configure owner input/select visibility and required attributes.
 */
function configureOwnerFields(ownerInput, ownerSelect, useAutocomplete, hasCategory) {
  if (ownerInput) {
    ownerInput.style.display = useAutocomplete ? "none" : ""
    ownerInput.required = !useAutocomplete
    ownerInput.disabled = !hasCategory
    ownerInput.placeholder = hasCategory ? "Name of dye, filter, etc." : "Select a category first"
  }
  if (ownerSelect) {
    ownerSelect.style.display = useAutocomplete ? "" : "none"
    ownerSelect.required = useAutocomplete
    ownerSelect.disabled = !hasCategory
  }
}

/**
 * Check if a select element has Select2 initialized.
 */
function isSelect2Initialized(selectEl) {
  return selectEl && window.$ && $(selectEl).hasClass("select2-hidden-accessible")
}

/**
 * Destroy Select2 instance if initialized.
 */
function destroySelect2(selectEl) {
  if (isSelect2Initialized(selectEl)) {
    $(selectEl).select2("destroy")
  }
}

/**
 * Initialize Select2 for protein autocomplete.
 */
function initOwnerSelect2(ownerSelect, category, initialValue, spectrum, el, state) {
  const url = AUTOCOMPLETE_URLS[category]
  if (!url || !window.$ || !ownerSelect) return

  destroySelect2(ownerSelect)

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

  // Pre-populate Select2 with initial value if provided
  if (initialValue) {
    const option = new Option(initialValue, initialValue, true, true)
    $(ownerSelect).append(option).trigger("change")
  }
}

/**
 * Update the owner warning based on current spectrum state.
 */
async function updateOwnerWarning(spectrum, index) {
  const warningEl = document.getElementById(`owner-warning-${index}`)
  if (!warningEl) return

  if (spectrum.owner && spectrum.category && spectrum.subtype) {
    spectrum.hasExactMatch = await checkSimilarOwners(
      spectrum.owner,
      spectrum.category,
      spectrum.subtype,
      warningEl
    )
  } else {
    spectrum.hasExactMatch = false
    warningEl.innerHTML = ""
    warningEl.classList.add("d-none")
  }
}

// ============================================================================
// Card Visibility Updates
// ============================================================================

function updateCardVisibility(index, { shouldNormalize, showFluor }) {
  const peakBadge = document.getElementById(`peak-badge-${index}`)
  const scaleContainer = document.getElementById(`scale-factor-container-${index}`)
  const optionalRow = document.querySelector(`.optional-fields-${index}`)
  const chartContainer = document.getElementById(`chart-container-${index}`)
  const chartHint = document.getElementById(`chart-hint-${index}`)

  if (peakBadge) peakBadge.style.display = shouldNormalize ? "" : "none"
  if (scaleContainer) scaleContainer.style.display = shouldNormalize ? "" : "none"
  if (optionalRow) optionalRow.style.display = shouldNormalize || showFluor ? "" : "none"
  if (chartContainer) chartContainer.style.cursor = shouldNormalize ? "crosshair" : "default"
  if (chartHint) chartHint.style.display = shouldNormalize ? "" : "none"

  document.querySelectorAll(`.fluor-field-${index}`).forEach((field) => {
    field.style.display = showFluor ? "" : "none"
  })
}

function updateScaleFactorUnits(index, subtype) {
  const unitsEl = document.getElementById(`scale-factor-units-${index}`)
  if (unitsEl) unitsEl.textContent = SCALE_UNITS[subtype] || ""
}

// ============================================================================
// Spectrum Processing
// ============================================================================

function processSpectrum(spectrum, index, el, state) {
  const { shouldNormalize } = getCategoryFlags(spectrum.category)
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
    // Use manualPeakWave if set (e.g., from restored state), otherwise use computed peak
    peakWave = spectrum.manualPeakWave ?? result.peakWave
  } else {
    // Use absolute values (no normalization)
    processedData = spectrum.interpolated
  }

  spectrum.processed = processedData

  // Update peak badge
  updatePeakBadge(index, shouldNormalize, peakWave)

  // Create or update chart
  updateSpectrumChart(spectrum, index, processedData, shouldNormalize, peakWave, el, state)
}

function updatePeakBadge(index, shouldNormalize, peakWave) {
  const peakBadge = document.getElementById(`peak-badge-${index}`)
  if (!peakBadge) return

  if (shouldNormalize) {
    peakBadge.style.display = ""
    peakBadge.textContent = peakWave !== null ? `Peak: ${peakWave} nm` : "No peak"
    peakBadge.className = `badge ${peakWave !== null ? "bg-primary" : "bg-secondary"}`
  } else {
    peakBadge.style.display = "none"
  }
}

function updateSpectrumChart(spectrum, index, processedData, shouldNormalize, peakWave, el, state) {
  const chartContainer = document.getElementById(`chart-container-${index}`)
  if (!chartContainer) return

  const chartOptions = {
    name: spectrum.columnName,
    normalized: shouldNormalize,
    rawData: shouldNormalize ? spectrum.interpolated : null,
    onClick: shouldNormalize
      ? (wavelength) => {
          spectrum.manualPeakWave = Math.round(wavelength)
          processSpectrum(spectrum, index, el, state)
          updateFormState(el, state)
        }
      : null,
  }

  // Track whether we had a click handler before
  const hadClickHandler = spectrum._hadClickHandler ?? false
  const hasClickHandler = chartOptions.onClick !== null
  spectrum._hadClickHandler = hasClickHandler

  // Recreate chart if click handler state changed, otherwise just update
  if (!spectrum.chartController) {
    spectrum.chartController = createSpectrumChart(chartContainer, processedData, chartOptions)
  } else if (hadClickHandler !== hasClickHandler) {
    // Click handler state changed - need to recreate chart
    spectrum.chartController.destroy()
    spectrum.chartController = createSpectrumChart(chartContainer, processedData, chartOptions)
  } else {
    spectrum.chartController.updateData(processedData, spectrum.columnName)
    spectrum.chartController.updateYAxis(
      processedData,
      shouldNormalize,
      shouldNormalize ? spectrum.interpolated : null
    )
  }

  // Update peak marker
  if (shouldNormalize && peakWave !== null) {
    spectrum.chartController.setPeakMarker(peakWave)
  } else {
    spectrum.chartController.clearAnnotations()
  }
}

// ============================================================================
// Form State Management
// ============================================================================

/**
 * Check for duplicate spectra within the form (same category, owner, subtype).
 * @param {Spectrum[]} spectra - Array of spectrum objects
 * @returns {Array<{indices: number[], key: string, owner: string, category: string, subtype: string}>}
 */
function checkForDuplicatesInForm(spectra) {
  const seen = new Map()
  const duplicates = []

  for (let i = 0; i < spectra.length; i++) {
    const s = spectra[i]
    if (!s.category || !s.subtype || !s.owner.trim()) continue

    const key = `${s.category}|${s.owner.trim().toLowerCase()}|${s.subtype}`

    if (seen.has(key)) {
      const firstIndex = seen.get(key)
      const existing = duplicates.find((d) => d.key === key)
      if (existing) {
        existing.indices.push(i)
      } else {
        duplicates.push({
          key,
          indices: [firstIndex, i],
          owner: s.owner,
          category: s.category,
          subtype: s.subtype,
        })
      }
    } else {
      seen.set(key, i)
    }
  }

  return duplicates
}

function updateFormState(el, state) {
  const validSpectra = state.spectra.filter(
    (s) => s.processed?.length > 0 && s.category && s.subtype && s.owner.trim()
  )

  /** @type {SpectrumJSON[]} */
  const spectraJson = validSpectra.map((s) => ({
    data: s.processed,
    category: s.category,
    owner: s.owner,
    subtype: s.subtype,
    scale_factor: s.scaleFactor,
    ph: FLUOR_CATEGORIES.has(s.category) ? s.ph : null,
    solvent: FLUOR_CATEGORIES.has(s.category) ? s.solvent : null,
    peak_wave: s.manualPeakWave ?? getPeakWave(s.processed),
    column_name: s.columnName,
  }))

  el.spectraJson.value = JSON.stringify(spectraJson)

  // Check for duplicates within this form submission
  const duplicatesInForm = checkForDuplicatesInForm(state.spectra)
  const duplicateIndices = new Set(duplicatesInForm.flatMap((d) => d.indices))

  // Track missing fields per spectrum
  const missing = { owners: [], categories: [], subtypes: [] }

  for (const [i, s] of state.spectra.entries()) {
    const hasOwner = s.owner.trim() !== ""
    const hasCategory = s.category !== ""
    const hasSubtype = s.subtype !== ""
    const isComplete = hasOwner && hasCategory && hasSubtype
    const isDuplicate = duplicateIndices.has(i)

    updateStatusIcon(i, isComplete, s.hasExactMatch, isDuplicate)
    updateFieldLabels(i, { hasOwner, hasCategory, hasSubtype })

    const spectrumLabel = `Spectrum ${i + 1}`
    if (!hasCategory) missing.categories.push(spectrumLabel)
    else if (!hasSubtype) missing.subtypes.push(spectrumLabel)
    if (!hasOwner) missing.owners.push(spectrumLabel)
  }

  // Check source/reference validation
  const hasSource = el.sourceInput?.value.trim() !== ""
  const hasReference = el.referenceInput?.value.trim() !== ""
  const hasValidReference = !hasReference || CONFIG.doiPattern.test(el.referenceInput.value.trim())
  const hasSourceOrRef = hasSource || hasReference

  // Check for exact matches in database (blocks submission)
  const hasAnyExactMatch = state.spectra.some((s) => s.hasExactMatch)

  // Check confirmation checkbox
  const isConfirmed = el.confirmationCheckbox?.checked ?? false

  // Update submit button
  const allComplete =
    missing.owners.length === 0 && missing.categories.length === 0 && missing.subtypes.length === 0
  const isValid =
    state.spectra.length > 0 &&
    allComplete &&
    hasSourceOrRef &&
    hasValidReference &&
    !hasAnyExactMatch &&
    duplicatesInForm.length === 0 &&
    isConfirmed

  el.submitBtn.disabled = !isValid

  updateValidationMessage(el, state, {
    missing,
    hasSourceOrRef,
    hasValidReference,
    hasAnyExactMatch,
    duplicatesInForm,
    isConfirmed,
  })
}

function updateStatusIcon(index, isComplete, hasExactMatch, isDuplicate) {
  const cardHeader = document.querySelector(`#spectrum-card-${index} .card-header`)
  let statusIcon = document.getElementById(`status-icon-${index}`)

  if (!statusIcon && cardHeader) {
    statusIcon = document.createElement("span")
    statusIcon.id = `status-icon-${index}`
    statusIcon.className = "me-2"
    cardHeader.querySelector("strong")?.prepend(statusIcon)
  }

  if (!statusIcon) return

  if (hasExactMatch) {
    statusIcon.innerHTML = '<span class="text-danger">✕</span> '
  } else if (isDuplicate) {
    statusIcon.innerHTML = '<span class="text-danger">⚠️</span> '
  } else {
    const icon = isComplete ? "✅" : "⚠️"
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

function updateValidationMessage(
  el,
  state,
  { missing, hasSourceOrRef, hasValidReference, hasAnyExactMatch, duplicatesInForm, isConfirmed }
) {
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

  if (hasAnyExactMatch) {
    const exactMatchNames = state.spectra
      .filter((s) => s.hasExactMatch)
      .map((s) => `"${escapeHtml(s.columnName)}"`)
      .join(", ")
    issues.push(
      `🚫 <strong>Exact match found:</strong> Cannot submit duplicate spectra: ${exactMatchNames}`
    )
  }

  if (duplicatesInForm.length > 0) {
    const dupMessages = duplicatesInForm.map((dup) => {
      const spectraNums = dup.indices.map((i) => i + 1).join(", ")
      return `Spectra ${spectraNums}: ${escapeHtml(dup.owner)} (${dup.subtype})`
    })
    issues.push(
      `🚫 <strong>Duplicate spectra in form:</strong> The following spectra have identical owner, category, and subtype:<br>${dupMessages.join("<br>")}`
    )
  }

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
      `<strong>${n}</strong> ${n === 1 ? "spectrum needs" : "spectra need"} an owner name: ${names}`
    )
  }

  if (!hasSourceOrRef) {
    issues.push(
      "Please provide at least one of <strong>Source</strong> or <strong>Primary Reference</strong>"
    )
  } else if (!hasValidReference) {
    issues.push("Primary Reference must be a valid DOI (e.g., 10.1234/example)")
  }

  if (!isConfirmed) {
    issues.push("Please check the confirmation box to verify data validity")
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
// State Restoration
// ============================================================================

/**
 * Restore UI state from spectra_json hidden field.
 * Called on page load to restore state after a form validation error.
 */
function restoreStateFromJson(el, state) {
  const jsonValue = el.spectraJson?.value?.trim()
  if (!jsonValue) return

  let spectraData
  try {
    spectraData = JSON.parse(jsonValue)
  } catch {
    return // Invalid JSON, nothing to restore
  }

  if (!Array.isArray(spectraData) || spectraData.length === 0) return

  // Rebuild UI from saved spectra data
  el.spectraPreview.style.display = "block"
  el.spectraPreview.innerHTML = ""
  el.spectraPreview.appendChild(createStep3Instructions())

  for (let i = 0; i < spectraData.length; i++) {
    const specData = spectraData[i]

    // Reconstruct spectrum object from saved data
    const spectrum = {
      columnName: specData.column_name || `Spectrum ${i + 1}`,
      raw: specData.data,
      interpolated: specData.data,
      processed: specData.data,
      chartController: null,
      category: specData.category || "",
      owner: specData.owner || "",
      subtype: specData.subtype || "",
      scaleFactor: specData.scale_factor,
      ph: specData.ph,
      solvent: specData.solvent || "",
      manualPeakWave: specData.peak_wave,
      hasExactMatch: false,
    }

    const card = createSpectrumCard(el, state, spectrum, i)
    el.spectraPreview.appendChild(card)
    processSpectrum(spectrum, i, el, state)
    state.spectra.push(spectrum)
  }

  el.globalSourceFields.style.display = "block"
  el.confirmationSection.style.display = "block"
  updateFormState(el, state)
}

// ============================================================================
// Utilities
// ============================================================================

function createStep3Instructions() {
  const div = document.createElement("div")
  div.className = "alert alert-info mb-3"
  div.innerHTML = `
    <strong>Step 3: Configure each spectrum</strong><br>
    <span class="text-primary"><i class="bi bi-1-circle me-1"></i>Select a <strong>Category</strong> and <strong>Subtype</strong> for each spectrum</span><br>
    <span class="text-info"><i class="bi bi-2-circle me-1"></i>Enter the <strong>Owner</strong> (protein, dye, filter name, etc.)</span><br>
    <span class="text-muted"><i class="bi bi-3-circle me-1"></i>Optionally click the chart to adjust peak location</span>
  `
  return div
}

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

/**
 * Form Controller Module
 *
 * Main orchestration for the spectrum submission form V2.
 * Manages the workflow from file upload through submission.
 */

import { renderColumnPicker } from "./column-picker.js"
import { extractSpectrum, parseCSV } from "./csv-parser.js"
import { interpolateToOneNm, normalize2P, normalizeSpectrum } from "./normalization.js"
import { createSpectrumChart } from "./spectrum-chart.js"

// Spectrum subtypes that need special 2P normalization
const TWO_PHOTON_SUBTYPES = ["2p"]

// Categories that should NOT be normalized (data is already percentage-based)
const NO_NORMALIZE_CATEGORIES = ["f", "c"] // Filter, Camera (lights ARE normalized)

// Categories that show pH/Solvent fields (biological spectra)
const BIO_CATEGORIES = ["d", "p"] // Dye, Protein

// Categories that use Select2 autocomplete for owner field (must select existing)
const AUTOCOMPLETE_CATEGORIES = {
  p: "/autocomplete-protein/", // Protein - must select existing
  // Filter, Dye, Camera, Light - allow free text
}

// Category to valid subtypes mapping (from backend Spectrum model)
const CATEGORY_SUBTYPES = {
  d: [
    // Dye
    { value: "ex", label: "Excitation" },
    { value: "ab", label: "Absorption" },
    { value: "em", label: "Emission" },
    { value: "2p", label: "Two-Photon" },
  ],
  p: [
    // Protein
    { value: "ex", label: "Excitation" },
    { value: "ab", label: "Absorption" },
    { value: "em", label: "Emission" },
    { value: "2p", label: "Two-Photon" },
  ],
  f: [
    // Filter
    { value: "bp", label: "Bandpass" },
    { value: "bx", label: "Bandpass-Ex" },
    { value: "bm", label: "Bandpass-Em" },
    { value: "sp", label: "Shortpass" },
    { value: "lp", label: "Longpass" },
    { value: "bs", label: "Beamsplitter" },
  ],
  c: [
    // Camera
    { value: "qe", label: "Quantum Efficiency" },
  ],
  l: [
    // Light source
    { value: "pd", label: "Power Distribution" },
  ],
}

// Category labels for display
const CATEGORY_LABELS = {
  d: "Dye",
  p: "Protein",
  f: "Filter",
  c: "Camera",
  l: "Light Source",
}

// Peak search radius in nm (±25nm around click)
const PEAK_SEARCH_RADIUS = 25

// DOI regex pattern
const DOI_PATTERN = /^10\.\d{4,}\/[^\s]+$/

// Scale factor units based on subtype
const SCALE_FACTOR_UNITS = {
  ex: "EC (M⁻¹cm⁻¹)",
  ab: "EC (M⁻¹cm⁻¹)",
  em: "QE (0-1)",
  "2p": "Cross Section (GM)",
  pd: "",
}

/**
 * Get scale factor units for a subtype.
 */
function getScaleFactorUnits(subtype) {
  return SCALE_FACTOR_UNITS[subtype] || ""
}

/**
 * Initialize the spectrum form.
 */
export function initSpectrumForm() {
  // Get DOM elements
  const form = document.getElementById("spectrum-form-v2")
  if (!form) {
    console.error("Spectrum form V2 not found")
    return
  }

  const fileInput = document.getElementById("id_file")
  const columnPickerContainer = document.getElementById("column-picker-container")
  const spectraPreviewContainer = document.getElementById("spectra-preview-container")
  const spectraJsonInput = document.getElementById("id_spectra_json")
  const submitBtn = form.querySelector('button[type="submit"]')
  const globalSourceFields = document.getElementById("global-source-fields")
  const confirmationSection = document.getElementById("confirmation-section")
  const sourceInput = document.getElementById("id_source")
  const referenceInput = document.getElementById("id_reference")

  // Global state
  const state = {
    parsedCSV: null,
    wavelengthCol: null,
    dataCols: [],
    spectra: [], // Array of spectrum states with all metadata
  }

  // DOI validation on input
  referenceInput?.addEventListener("input", (e) => {
    const value = e.target.value.trim()
    if (value && !DOI_PATTERN.test(value)) {
      e.target.setCustomValidity("Please enter a valid DOI (e.g., 10.1234/example)")
    } else {
      e.target.setCustomValidity("")
    }
    updateFormState()
  })

  // Source input - trigger validation
  sourceInput?.addEventListener("input", () => {
    updateFormState()
  })

  // File upload handler
  fileInput?.addEventListener("change", async (e) => {
    const file = e.target.files[0]
    if (!file) return

    try {
      const text = await file.text()
      state.parsedCSV = parseCSV(text)

      // Show column picker
      renderColumnPicker(columnPickerContainer, state.parsedCSV, (waveCol, dataCols) => {
        state.wavelengthCol = waveCol
        state.dataCols = dataCols
        processSelectedColumns()
      })

      // Hide other sections until columns are selected
      spectraPreviewContainer.style.display = "none"
      spectraPreviewContainer.innerHTML = ""
      globalSourceFields.style.display = "none"
      confirmationSection.style.display = "none"

      // Disable submit until processed
      submitBtn.disabled = true
    } catch (error) {
      console.error("Error parsing file:", error)
      showAlert(columnPickerContainer, `Error parsing file: ${error.message}`, "danger")
    }
  })

  /**
   * Process selected columns and create spectrum previews.
   */
  function processSelectedColumns() {
    // Clear previous spectra
    state.spectra.forEach((s) => {
      s.chartController?.destroy()
    })
    state.spectra = []

    // Hide column picker, show preview container
    columnPickerContainer.style.display = "none"
    spectraPreviewContainer.style.display = "block"
    spectraPreviewContainer.innerHTML = ""

    // Process each selected data column
    state.dataCols.forEach((dataColIndex, i) => {
      const columnName = state.parsedCSV.headers[dataColIndex]
      const rawData = extractSpectrum(state.parsedCSV, state.wavelengthCol, dataColIndex)

      if (rawData.length < 2) {
        showAlert(
          spectraPreviewContainer,
          `Column "${columnName}" has insufficient data points`,
          "warning"
        )
        return
      }

      // Interpolate once and store
      const interpolated = interpolateToOneNm(rawData)

      // Create spectrum entry with all metadata
      const spectrumState = {
        columnName,
        raw: rawData,
        interpolated,
        processed: null,
        chartController: null,
        // Per-spectrum metadata
        category: "", // Must be selected by user
        owner: "",
        subtype: "", // Set when category is selected
        scaleFactor: null,
        ph: null,
        solvent: "",
        manualPeakWave: null,
      }

      // Create preview card
      const card = createSpectrumCard(spectrumState, i)
      spectraPreviewContainer.appendChild(card)

      // Initial processing
      processSpectrum(spectrumState, i)

      state.spectra.push(spectrumState)
    })

    // Show global source fields and confirmation
    globalSourceFields.style.display = "block"
    confirmationSection.style.display = "block"

    // Update form state
    updateFormState()
  }

  /**
   * Create a preview card for a spectrum with all metadata fields.
   */
  function createSpectrumCard(spectrumState, index) {
    const card = document.createElement("div")
    card.className = "card mb-3 spectrum-card"
    card.id = `spectrum-card-${index}`

    const minWave = Math.round(spectrumState.raw[0][0])
    const maxWave = Math.round(spectrumState.raw[spectrumState.raw.length - 1][0])

    // Build category options with placeholder
    const categoryOptions =
      `<option value="" ${!spectrumState.category ? "selected" : ""}>-----</option>` +
      Object.entries(CATEGORY_LABELS)
        .map(
          ([value, label]) =>
            `<option value="${value}" ${value === spectrumState.category ? "selected" : ""}>${label}</option>`
        )
        .join("")

    // Build subtype options
    const subtypeOptions = buildSubtypeOptions(spectrumState.category, spectrumState.subtype)

    // Check if should normalize and show bio fields
    const shouldNormalize = !NO_NORMALIZE_CATEGORIES.includes(spectrumState.category)
    const showBioFields = BIO_CATEGORIES.includes(spectrumState.category)
    const useAutocomplete = spectrumState.category in AUTOCOMPLETE_CATEGORIES

    card.innerHTML = `
      <div class="card-header d-flex justify-content-between align-items-center">
        <strong>Spectrum ${index + 1}: ${escapeHtml(spectrumState.columnName)}</strong>
        <span class="badge bg-secondary" id="peak-badge-${index}" style="${shouldNormalize ? "" : "display: none;"}">Peak: --</span>
      </div>
      <div class="card-body">
        <!-- Row 1: Category and Subtype -->
        <div class="row mb-3">
          <div class="col-md-6">
            <label class="form-label" id="category-label-${index}">Spectrum Category <span class="text-danger">*</span></label>
            <select class="form-select form-select-sm" id="category-select-${index}" required>
              ${categoryOptions}
            </select>
          </div>
          <div class="col-md-6">
            <label class="form-label" id="subtype-label-${index}">Subtype <span class="text-danger">*</span></label>
            <select class="form-select form-select-sm" id="subtype-select-${index}" required>
              ${subtypeOptions}
            </select>
          </div>
        </div>

        <!-- Row 2: Owner (full width) -->
        <div class="row mb-3">
          <div class="col-12">
            <label class="form-label" id="owner-label-${index}">Owner <span class="text-danger">*</span></label>
            <!-- Text input for free-text categories -->
            <input type="text" class="form-control form-control-sm" id="owner-input-${index}"
                   placeholder="Name of dye, filter, etc." required
                   style="${useAutocomplete ? "display: none;" : ""}">
            <!-- Select for autocomplete categories (Protein) -->
            <select class="form-select form-select-sm" id="owner-select-${index}"
                    style="${useAutocomplete ? "" : "display: none;"}">
              <option value="">Select a protein...</option>
            </select>
          </div>
        </div>

        <!-- Row 3: Scale Factor, pH, Solvent (conditional) -->
        <div class="row mb-3 optional-fields-${index}" style="${shouldNormalize || showBioFields ? "" : "display: none;"}">
          <div class="col-md-4" id="scale-factor-container-${index}" style="${shouldNormalize ? "" : "display: none;"}">
            <label class="form-label">Scale Factor <small class="text-muted">(optional)</small></label>
            <div class="d-flex align-items-center gap-2">
              <input type="number" class="form-control form-control-sm no-spinners flex-grow-1" id="scale-factor-${index}" step="any">
              <small class="text-muted text-nowrap flex-shrink-0" id="scale-factor-units-${index}">${getScaleFactorUnits(spectrumState.subtype)}</small>
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

        <!-- Chart -->
        <div id="chart-container-${index}" style="cursor: ${shouldNormalize ? "crosshair" : "default"};"></div>

        <div class="small text-muted" id="chart-hint-${index}" style="${shouldNormalize ? "" : "display: none;"}">
          <em>Click on the chart to set peak location (searches ±${PEAK_SEARCH_RADIUS}nm around click)</em>
        </div>

        <!-- Footer -->
        <div class="mt-2 d-flex justify-content-between align-items-center">
          <span class="text-muted small">
            ${spectrumState.raw.length} points | ${minWave}-${maxWave} nm
          </span>
          <button type="button" class="btn btn-sm btn-outline-danger" id="remove-btn-${index}">
            Remove
          </button>
        </div>
      </div>
    `

    // Add event listeners after DOM insertion
    setTimeout(() => {
      const categorySelect = document.getElementById(`category-select-${index}`)
      const ownerInput = document.getElementById(`owner-input-${index}`)
      const ownerSelect = document.getElementById(`owner-select-${index}`)
      const subtypeSelect = document.getElementById(`subtype-select-${index}`)
      const scaleFactorInput = document.getElementById(`scale-factor-${index}`)
      const phInput = document.getElementById(`ph-input-${index}`)
      const solventInput = document.getElementById(`solvent-input-${index}`)
      const removeBtn = document.getElementById(`remove-btn-${index}`)

      // Initialize Select2 for autocomplete categories
      function initOwnerSelect2(category) {
        const url = AUTOCOMPLETE_CATEGORIES[category]
        if (!url || !window.$ || !ownerSelect) return

        // Destroy existing Select2 if any
        if ($(ownerSelect).hasClass("select2-hidden-accessible")) {
          $(ownerSelect).select2("destroy")
        }

        $(ownerSelect).select2({
          theme: "bootstrap-5",
          placeholder: "Search for a protein...",
          allowClear: true,
          ajax: {
            url: url,
            dataType: "json",
            delay: 250,
            cache: true,
            data: (params) => ({ q: params.term, page: params.page }),
            processResults: (data) => data,
          },
        })

        // Handle selection change
        $(ownerSelect).on("select2:select select2:clear", () => {
          const selected = $(ownerSelect).select2("data")
          spectrumState.owner = selected?.[0]?.text || ""
          updateFormState()
        })
      }

      // Initialize if starting with autocomplete category
      if (spectrumState.category in AUTOCOMPLETE_CATEGORIES) {
        initOwnerSelect2(spectrumState.category)
      }

      // Category change - update subtypes and show/hide fields
      categorySelect?.addEventListener("change", (e) => {
        const newCategory = e.target.value
        spectrumState.category = newCategory

        // Update subtype options (auto-select if only one option, else reset to placeholder)
        const subtypes = CATEGORY_SUBTYPES[newCategory] || []
        spectrumState.subtype = subtypes.length === 1 ? subtypes[0].value : ""
        subtypeSelect.innerHTML = buildSubtypeOptions(newCategory, spectrumState.subtype)

        const shouldNorm = !NO_NORMALIZE_CATEGORIES.includes(newCategory)
        const showBio = BIO_CATEGORIES.includes(newCategory)
        const useAutocomplete = newCategory in AUTOCOMPLETE_CATEGORIES

        // Switch between text input and select for owner
        if (ownerInput) ownerInput.style.display = useAutocomplete ? "none" : ""
        if (ownerSelect) ownerSelect.style.display = useAutocomplete ? "" : "none"

        // Reset owner value when switching modes
        spectrumState.owner = ""
        if (ownerInput) ownerInput.value = ""
        if (ownerSelect && $(ownerSelect).hasClass("select2-hidden-accessible")) {
          $(ownerSelect).val(null).trigger("change")
        }

        // Initialize or destroy Select2
        if (useAutocomplete) {
          initOwnerSelect2(newCategory)
        } else if (
          ownerSelect &&
          window.$ &&
          $(ownerSelect).hasClass("select2-hidden-accessible")
        ) {
          $(ownerSelect).select2("destroy")
        }

        // Show/hide peak badge
        const peakBadge = document.getElementById(`peak-badge-${index}`)
        if (peakBadge) peakBadge.style.display = shouldNorm ? "" : "none"

        // Show/hide bio fields (pH, Solvent)
        document.querySelectorAll(`.bio-field-${index}`).forEach((el) => {
          el.style.display = showBio ? "" : "none"
        })

        // Show/hide scale factor and update units
        const scaleContainer = document.getElementById(`scale-factor-container-${index}`)
        if (scaleContainer) scaleContainer.style.display = shouldNorm ? "" : "none"
        const unitsEl = document.getElementById(`scale-factor-units-${index}`)
        if (unitsEl) unitsEl.textContent = getScaleFactorUnits(spectrumState.subtype)

        // Show/hide the entire optional fields row
        const optionalRow = document.querySelector(`.optional-fields-${index}`)
        if (optionalRow) {
          optionalRow.style.display = shouldNorm || showBio ? "" : "none"
        }

        // Update chart cursor and hint
        const chartContainer = document.getElementById(`chart-container-${index}`)
        const chartHint = document.getElementById(`chart-hint-${index}`)
        if (chartContainer) chartContainer.style.cursor = shouldNorm ? "crosshair" : "default"
        if (chartHint) chartHint.style.display = shouldNorm ? "" : "none"

        // Reset manual peak
        spectrumState.manualPeakWave = null

        // Reprocess
        processSpectrum(spectrumState, index)
        updateFormState()
      })

      ownerInput?.addEventListener("input", (e) => {
        spectrumState.owner = e.target.value.trim()
        updateFormState()
      })

      subtypeSelect?.addEventListener("change", (e) => {
        spectrumState.subtype = e.target.value
        spectrumState.manualPeakWave = null

        // Update scale factor units
        const unitsEl = document.getElementById(`scale-factor-units-${index}`)
        if (unitsEl) unitsEl.textContent = getScaleFactorUnits(spectrumState.subtype)

        processSpectrum(spectrumState, index)
        updateFormState()
      })

      scaleFactorInput?.addEventListener("change", (e) => {
        const val = parseFloat(e.target.value)
        spectrumState.scaleFactor = Number.isNaN(val) ? null : val
        updateFormState()
      })

      phInput?.addEventListener("change", (e) => {
        const val = parseFloat(e.target.value)
        spectrumState.ph = Number.isNaN(val) ? null : val
        updateFormState()
      })

      solventInput?.addEventListener("input", (e) => {
        spectrumState.solvent = e.target.value.trim()
        updateFormState()
      })

      removeBtn?.addEventListener("click", () => {
        // Clean up Select2 if initialized
        if (ownerSelect && window.$ && $(ownerSelect).hasClass("select2-hidden-accessible")) {
          $(ownerSelect).select2("destroy")
        }
        spectrumState.chartController?.destroy()
        card.remove()
        state.spectra = state.spectra.filter((s) => s !== spectrumState)
        updateFormState()

        // Hide global fields if no spectra left
        if (state.spectra.length === 0) {
          globalSourceFields.style.display = "none"
          confirmationSection.style.display = "none"
        }
      })
    }, 0)

    return card
  }

  /**
   * Build subtype <option> HTML for a category.
   */
  function buildSubtypeOptions(category, selectedValue) {
    if (!category) {
      return `<option value="" selected>Select category first...</option>`
    }
    const subtypes = CATEGORY_SUBTYPES[category] || []
    // Skip placeholder if only one subtype option
    if (subtypes.length === 1) {
      return `<option value="${subtypes[0].value}" selected>${subtypes[0].label}</option>`
    }
    const placeholder = `<option value="" ${!selectedValue ? "selected" : ""}>-----</option>`
    return (
      placeholder +
      subtypes
        .map(
          ({ value, label }) =>
            `<option value="${value}" ${value === selectedValue ? "selected" : ""}>${label}</option>`
        )
        .join("")
    )
  }

  /**
   * Process a spectrum (normalize, update chart).
   */
  function processSpectrum(spectrumState, index) {
    const interpolated = spectrumState.interpolated
    const shouldNormalize = !NO_NORMALIZE_CATEGORIES.includes(spectrumState.category)

    let processedData
    let peakWave = null

    if (shouldNormalize) {
      // Normalize based on subtype and manual peak
      let result
      if (spectrumState.manualPeakWave != null) {
        const rangeMin = spectrumState.manualPeakWave - PEAK_SEARCH_RADIUS
        const rangeMax = spectrumState.manualPeakWave + PEAK_SEARCH_RADIUS
        const options = { rangeMin, rangeMax }

        if (TWO_PHOTON_SUBTYPES.includes(spectrumState.subtype)) {
          result = normalize2P(interpolated, options)
        } else {
          result = normalizeSpectrum(interpolated, options)
        }
      } else {
        if (TWO_PHOTON_SUBTYPES.includes(spectrumState.subtype)) {
          result = normalize2P(interpolated, {})
        } else {
          result = normalizeSpectrum(interpolated, {})
        }
      }
      processedData = result.normalized
      peakWave = result.peakWave
    } else {
      // For filters/cameras: use data as-is
      const maxVal = Math.max(...interpolated.map(([_, y]) => y))
      if (maxVal > 1) {
        processedData = interpolated.map(([x, y]) => [x, y / 100])
      } else {
        processedData = interpolated
      }
    }

    spectrumState.processed = processedData

    // Update peak badge
    const peakBadge = document.getElementById(`peak-badge-${index}`)
    if (peakBadge) {
      if (peakWave !== null) {
        peakBadge.textContent = `Peak: ${peakWave} nm`
        peakBadge.className = "badge bg-primary"
      } else {
        peakBadge.textContent = "No peak"
        peakBadge.className = "badge bg-secondary"
      }
    }

    // Create or update chart
    const chartContainer = document.getElementById(`chart-container-${index}`)
    const chartOptions = {
      name: spectrumState.columnName,
      onClick: shouldNormalize
        ? (wavelength) => handleChartClick(spectrumState, index, wavelength)
        : null,
    }

    if (!spectrumState.chartController && chartContainer) {
      spectrumState.chartController = createSpectrumChart(
        chartContainer,
        processedData,
        chartOptions
      )
      if (peakWave !== null) {
        spectrumState.chartController.setPeakMarker(peakWave)
      }
    } else if (spectrumState.chartController) {
      spectrumState.chartController.updateData(processedData, spectrumState.columnName)
      if (peakWave !== null) {
        spectrumState.chartController.setPeakMarker(peakWave)
      } else {
        spectrumState.chartController.clearAnnotations()
      }
    }
  }

  /**
   * Handle click on chart to set manual peak location.
   */
  function handleChartClick(spectrumState, index, wavelength) {
    spectrumState.manualPeakWave = Math.round(wavelength)
    processSpectrum(spectrumState, index)
    updateFormState()
  }

  /**
   * Update the hidden JSON field with current spectra data.
   */
  function updateFormState() {
    // Only include spectra with required fields filled
    const validSpectra = state.spectra.filter(
      (s) =>
        s.processed && s.processed.length > 0 && s.category && s.subtype && s.owner.trim() !== ""
    )

    const spectraJson = validSpectra.map((s) => ({
      data: s.processed,
      category: s.category,
      owner: s.owner,
      subtype: s.subtype,
      scale_factor: s.scaleFactor,
      ph: BIO_CATEGORIES.includes(s.category) ? s.ph : null,
      solvent: BIO_CATEGORIES.includes(s.category) ? s.solvent : null,
      peak_wave: getPeakWave(s.processed),
      column_name: s.columnName,
    }))

    spectraJsonInput.value = JSON.stringify(spectraJson)

    // Update validation UI for each spectrum
    const missingFields = { owners: [], categories: [], subtypes: [] }
    state.spectra.forEach((s, index) => {
      const hasOwner = s.owner.trim() !== ""
      const hasCategory = s.category !== ""
      const hasSubtype = s.subtype !== ""
      const isComplete = hasOwner && hasCategory && hasSubtype
      const cardHeader = document.querySelector(`#spectrum-card-${index} .card-header`)
      const ownerLabel = document.getElementById(`owner-label-${index}`)
      const categoryLabel = document.getElementById(`category-label-${index}`)
      const subtypeLabel = document.getElementById(`subtype-label-${index}`)

      // Update card status indicator
      let statusIcon = document.getElementById(`status-icon-${index}`)
      if (!statusIcon && cardHeader) {
        statusIcon = document.createElement("span")
        statusIcon.id = `status-icon-${index}`
        statusIcon.className = "me-2"
        cardHeader.querySelector("strong").prepend(statusIcon)
      }
      if (statusIcon) {
        if (isComplete) {
          statusIcon.innerHTML = `<span class="text-success">✓</span> `
        } else {
          statusIcon.innerHTML = `<span class="text-warning">!</span> `
          if (!hasCategory) missingFields.categories.push(s.columnName)
          else if (!hasSubtype) missingFields.subtypes.push(s.columnName)
          if (!hasOwner) missingFields.owners.push(s.columnName)
        }
      }

      // Bold/unbold labels for incomplete fields
      if (ownerLabel) ownerLabel.style.fontWeight = hasOwner ? "normal" : "bold"
      if (categoryLabel) categoryLabel.style.fontWeight = hasCategory ? "normal" : "bold"
      if (subtypeLabel) subtypeLabel.style.fontWeight = hasSubtype ? "normal" : "bold"
    })

    // Check if at least one of source/reference is provided
    const hasSource = sourceInput?.value.trim() !== ""
    const hasReference = referenceInput?.value.trim() !== ""
    const hasValidReference =
      !referenceInput?.value.trim() || DOI_PATTERN.test(referenceInput.value.trim())
    const hasSourceOrRef = hasSource || hasReference

    // Update submit button and message
    const allComplete =
      missingFields.owners.length === 0 &&
      missingFields.categories.length === 0 &&
      missingFields.subtypes.length === 0
    const isValid = state.spectra.length > 0 && allComplete && hasSourceOrRef && hasValidReference
    submitBtn.disabled = !isValid

    // Show/update validation message
    updateValidationMessage(missingFields, hasSourceOrRef, hasValidReference)
  }

  /**
   * Update the validation message near the submit button.
   */
  function updateValidationMessage(missingFields, hasSourceOrRef, hasValidReference) {
    let messageEl = document.getElementById("validation-message")

    if (!messageEl) {
      messageEl = document.createElement("div")
      messageEl.id = "validation-message"
      messageEl.className = "mb-3"
      // Insert before the button container, not inside it
      const buttonContainer = submitBtn.parentElement
      buttonContainer.parentElement.insertBefore(messageEl, buttonContainer)
    }

    if (state.spectra.length === 0) {
      messageEl.innerHTML = ""
      return
    }

    // Collect all issues
    const issues = []

    if (missingFields.categories.length > 0) {
      const count = missingFields.categories.length
      const word = count === 1 ? "spectrum needs" : "spectra need"
      issues.push(
        `<strong>${count}</strong> ${word} a category: ${missingFields.categories.map((n) => `"${escapeHtml(n)}"`).join(", ")}`
      )
    }

    if (missingFields.subtypes.length > 0) {
      const count = missingFields.subtypes.length
      const word = count === 1 ? "spectrum needs" : "spectra need"
      issues.push(
        `<strong>${count}</strong> ${word} a subtype: ${missingFields.subtypes.map((n) => `"${escapeHtml(n)}"`).join(", ")}`
      )
    }

    if (missingFields.owners.length > 0) {
      const count = missingFields.owners.length
      const word = count === 1 ? "spectrum needs" : "spectra need"
      issues.push(
        `<strong>${count}</strong> ${word} an owner: ${missingFields.owners.map((n) => `"${escapeHtml(n)}"`).join(", ")}`
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
          <ul class="mb-0 ps-3">
            ${issues.map((i) => `<li>${i}</li>`).join("")}
          </ul>
        </div>
      `
    } else {
      messageEl.innerHTML = `
        <div class="alert alert-success py-2 mb-2">
          ✓ All spectra ready to submit
        </div>
      `
    }
  }

  /**
   * Get peak wavelength from normalized data.
   */
  function getPeakWave(data) {
    if (!data || data.length === 0) return null

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

  /**
   * Show an alert in a container.
   */
  function showAlert(container, message, type = "info") {
    const alert = document.createElement("div")
    alert.className = `alert alert-${type} alert-dismissible fade show`
    alert.innerHTML = `
      ${message}
      <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
    `
    container.prepend(alert)
  }

  /**
   * Escape HTML to prevent XSS.
   */
  function escapeHtml(str) {
    const div = document.createElement("div")
    div.textContent = str
    return div.innerHTML
  }
}

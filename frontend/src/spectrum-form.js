/**
 * Modern spectrum submission form with Alpine.js and Tom-Select
 *
 * Alpine.js handles reactive state management and form logic
 * Tom-Select provides the autocomplete widget
 */

import "tom-select/dist/css/tom-select.bootstrap4.css"
import Alpine from "@alpinejs/csp"
import { createAutocompleteSelect } from "./lib/autocomplete-select.js"

// Configuration for valid subtypes per category
const VALID_SUBTYPES = {
  d: ["ex", "ab", "em", "2p"], // Dye
  p: ["ex", "ab", "em", "2p"], // Protein
  l: ["pd"], // Light
  f: ["bp", "bx", "bm", "sp", "lp", "bs"], // Filter
  c: ["qe"], // Camera
  "": [],
}

// Register Alpine component
Alpine.data("spectrumForm", () => ({
  // State
  category: "",
  subtype: "",
  previewData: null,
  submitting: false,
  tomSelectInstance: null,

  // Initialize
  init() {
    // Get select elements by ID
    const categorySelect = document.getElementById("id_category")
    const subtypeSelect = document.getElementById("id_subtype")

    // Store original subtype options
    this.originalSubtypes = {}
    if (subtypeSelect) {
      Array.from(subtypeSelect.options).forEach((opt) => {
        if (opt.value) {
          this.originalSubtypes[opt.value] = opt.textContent
        }
      })
    }

    // Watch category changes to update UI
    this.$watch("category", (value) => {
      this.updateSubtypeOptions(value)
      this.updateOwnerFields(value)
    })

    // Set up two-way binding for category select
    if (categorySelect) {
      this.category = categorySelect.value
      categorySelect.addEventListener("change", (e) => {
        this.category = e.target.value
      })
    }

    // Set up two-way binding for subtype select
    if (subtypeSelect) {
      this.subtype = subtypeSelect.value
      subtypeSelect.addEventListener("change", (e) => {
        this.subtype = e.target.value
      })
    }
  },

  // Computed properties as getters
  get showProteinOwner() {
    return this.category === "p"
  },

  get showOtherOwner() {
    return this.category !== "p"
  },

  get showBioFields() {
    return this.category === "d" || this.category === "p"
  },

  get ownerTypeLabel() {
    const select = document.getElementById("id_category")
    return select?.options[select.selectedIndex]?.text || "Owner"
  },

  get hasPreview() {
    return this.previewData !== null
  },

  get submitButtonText() {
    if (this.submitting) return "Processing..."
    const activeTab = document.querySelector("#data-source-tabs .nav-link.active")
    const isFileTab = activeTab?.id === "file-tab"
    const fileInput = document.getElementById("id_file")
    const dataInput = document.getElementById("id_spectral_data")
    const hasData = isFileTab ? fileInput?.value : dataInput?.value?.trim()
    return hasData ? "Preview Spectrum" : "Submit"
  },

  // Methods
  updateSubtypeOptions(category) {
    const validSubtypes = VALID_SUBTYPES[category] || []
    const subtypeSelect = document.getElementById("id_subtype")

    if (!subtypeSelect) return

    // Remove all options except empty
    Array.from(subtypeSelect.options).forEach((opt) => {
      if (opt.value !== "") {
        opt.remove()
      }
    })

    // Add back only valid options
    validSubtypes.forEach((value) => {
      if (this.originalSubtypes[value]) {
        const option = document.createElement("option")
        option.value = value
        option.textContent = this.originalSubtypes[value]
        subtypeSelect.appendChild(option)
      }
    })

    // Reset subtype if current value is not valid
    if (!validSubtypes.includes(this.subtype)) {
      this.subtype = ""
    }
  },

  updateOwnerFields(category) {
    if (category === "p") {
      // Initialize Tom-Select for protein autocomplete
      this.initTomSelect()
    }
  },

  initTomSelect() {
    const select = this.$refs.ownerStateSelect
    const searchUrl = this.$el.dataset.stateSearchUrl

    if (!select || !searchUrl || this.tomSelectInstance) return

    this.tomSelectInstance = createAutocompleteSelect(select, {
      searchUrl,
      searchFields: ["text", "protein_name", "state_name"],
      placeholder: "Type to search proteins...",
      noResultsText: "No proteins found",
    })
  },

  async handleSubmit(event) {
    event.preventDefault()

    if (this.hasPreview) {
      // Already previewed, submit the form natively
      this.$refs.form.submit()
      return
    }

    // Show preview first
    await this.showPreview()
  },

  async showPreview() {
    this.submitting = true

    try {
      const formData = new FormData(this.$refs.form)

      // Determine active tab and set data source
      const activeTab = document.querySelector("#data-source-tabs .nav-link.active")
      const isManualTab = activeTab?.id === "manual-tab"

      if (isManualTab) {
        formData.delete("file")
        formData.append("data_source", "manual")
      } else {
        formData.set("spectral_data", "")
        formData.append("data_source", "file")
      }

      const previewUrl = this.$el.dataset.previewUrl

      const response = await fetch(previewUrl, {
        method: "POST",
        headers: {
          "X-CSRFToken": this.getCsrfToken(),
          "X-Requested-With": "XMLHttpRequest",
        },
        body: formData,
      })

      const data = await response.json()

      if (response.ok && data.success) {
        this.previewData = data.preview // Set this FIRST to trigger x-show

        // Use $nextTick to allow Alpine's reactivity to fully update the DOM
        this.$nextTick(() => {
          this.displayPreview(data)
        })
      } else {
        this.showError(data.error || "Failed to generate preview", data.details, data.form_errors)
      }
    } catch (error) {
      this.showError(`An unexpected error occurred: ${error.message}`)
    } finally {
      this.submitting = false
    }
  },

  displayPreview(response) {
    // Clear any previous errors
    this.$el.querySelectorAll(".alert-danger").forEach((alert) => {
      alert.remove()
    })

    const preview = response.preview

    // Add data_source field to form for final submission
    const activeTab = document.querySelector("#data-source-tabs .nav-link.active")
    const dataSource = activeTab?.id === "manual-tab" ? "manual" : "file"

    // Remove existing data_source hidden input and add new one
    this.$refs.form.querySelector('input[name="data_source"]')?.remove()
    const hiddenInput = Object.assign(document.createElement("input"), {
      type: "hidden",
      name: "data_source",
      value: dataSource,
    })
    this.$refs.form.appendChild(hiddenInput)

    // Helper to query within preview section (refs not available when x-show is false)
    const $ = (ref) => document.querySelector(`#spectrum-preview-section [x-ref='${ref}']`)

    // Update preview content
    const elements = {
      previewMessage: $("previewMessage"),
      previewPeakWave: $("previewPeakWave"),
      previewWaveRange: $("previewWaveRange"),
      previewDataPoints: $("previewDataPoints"),
      previewChart: $("previewChart"),
    }

    if (elements.previewMessage) elements.previewMessage.textContent = response.message
    if (elements.previewPeakWave) elements.previewPeakWave.textContent = preview.peak_wave || "N/A"
    if (elements.previewWaveRange)
      elements.previewWaveRange.textContent = preview.wave_range || "N/A"
    if (elements.previewDataPoints)
      elements.previewDataPoints.textContent = preview.data_points || "N/A"
    if (elements.previewChart) elements.previewChart.innerHTML = preview.svg || ""

    // Scroll to preview
    document.querySelector("#spectrum-preview-section")?.scrollIntoView({ behavior: "smooth" })
  },

  hidePreview() {
    this.previewData = null
  },

  showError(message, details, formErrors) {
    let errorHtml = '<div class="alert alert-danger alert-dismissible fade show" role="alert">'
    errorHtml += `<strong>Error:</strong> ${message}`

    if (details) {
      errorHtml += `<br><small>${details}</small>`
    }

    if (formErrors && Object.keys(formErrors).length > 0) {
      errorHtml += "<br><br><strong>Form Issues:</strong><ul>"
      for (const [field, errors] of Object.entries(formErrors)) {
        errorHtml += `<li><strong>${field}:</strong> ${errors.join(", ")}</li>`
      }
      errorHtml += "</ul>"
    }

    errorHtml += '<button type="button" class="close" data-dismiss="alert" aria-label="Close">'
    errorHtml += '<span aria-hidden="true">&times;</span></button></div>'

    this.$el.insertAdjacentHTML("afterbegin", errorHtml)
    this.$el.scrollIntoView({ behavior: "smooth", block: "start" })
  },

  getCsrfToken() {
    const name = "csrftoken"
    const cookies = document.cookie.split(";")
    for (const cookie of cookies) {
      const trimmed = cookie.trim()
      if (trimmed.startsWith(`${name}=`)) {
        return decodeURIComponent(trimmed.substring(name.length + 1))
      }
    }
    return null
  },
}))

// Start Alpine when DOM is ready
if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", () => {
    Alpine.start()
  })
} else {
  Alpine.start()
}

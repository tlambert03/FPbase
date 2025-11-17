/**
 * Progressive enhancement for modern spectrum submission form
 *
 * This module provides client-side enhancements without requiring
 * specific JavaScript libraries like jQuery, Select2, or autocomplete-light.
 * It works with plain vanilla JS and Bootstrap 4.
 */

(function () {
  "use strict";

  // State
  let currentPreviewData = null;

  // Configuration for valid subtypes per category
  const VALID_SUBTYPES = {
    d: ["ex", "ab", "em", "2p"], // Dye
    p: ["ex", "ab", "em", "2p"], // Protein
    l: ["pd"], // Light
    f: ["bp", "bx", "bm", "sp", "lp", "bs"], // Filter
    c: ["qe"], // Camera
    "": [],
  };

  // DOM elements (will be initialized on DOMContentLoaded)
  let elements = {};

  /**
   * Get CSRF token from cookie
   */
  function getCsrfToken() {
    const name = "csrftoken";
    const cookies = document.cookie.split(";");
    for (const cookie of cookies) {
      const trimmed = cookie.trim();
      if (trimmed.startsWith(name + "=")) {
        return decodeURIComponent(trimmed.substring(name.length + 1));
      }
    }
    return null;
  }

  /**
   * Initialize DOM element references
   */
  function initElements() {
    elements = {
      form: document.getElementById("spectrum-form"),
      categoryField: document.getElementById("id_category"),
      subtypeField: document.getElementById("id_subtype"),
      proteinOwnerField: document.getElementById("protein-owner-field"),
      otherOwnerField: document.getElementById("other-owner-field"),
      ownerStateSelect: document.getElementById("id_owner_state"),
      ownerInput: document.getElementById("id_owner"),
      ownerTypeLabel: document.getElementById("owner-type-label"),
      bioFieldsRow: document.getElementById("bio-fields-row"),
      phField: document.getElementById("id_ph"),
      solventField: document.getElementById("id_solvent"),
      fileInput: document.getElementById("id_file"),
      dataInput: document.getElementById("id_spectral_data"),
      submitBtn: document.getElementById("submit-btn"),
      clearStateBtn: document.getElementById("clear-state-btn"),
      formSection: document.getElementById("spectrum-form-section"),
      previewSection: document.getElementById("spectrum-preview-section"),
      previewMessage: document.getElementById("spectrum-preview-message"),
      previewChart: document.getElementById("spectrum-preview-chart"),
      previewPeakWave: document.getElementById("preview-peak-wave"),
      previewWaveRange: document.getElementById("preview-wave-range"),
      previewDataPoints: document.getElementById("preview-data-points"),
      editDataBtn: document.getElementById("edit-data-btn"),
      submitFinalBtn: document.getElementById("submit-final-btn"),
      fileTab: document.getElementById("file-tab"),
      manualTab: document.getElementById("manual-tab"),
    };
  }

  /**
   * Store original subtype options
   */
  let originalSubtypeOptions = null;
  function storeOriginalSubtypes() {
    if (originalSubtypeOptions) return;

    originalSubtypeOptions = {};
    const options = elements.subtypeField.querySelectorAll("option");
    options.forEach((opt) => {
      if (opt.value) {
        originalSubtypeOptions[opt.value] = opt.textContent;
      }
    });
  }

  /**
   * Handle category selection changes
   */
  function handleCategoryChange() {
    const category = elements.categoryField.value;

    // Update subtype options
    updateSubtypeOptions(category);

    // Update owner field visibility and labels
    updateOwnerFields(category);

    // Show/hide bio fields (pH, solvent)
    updateBioFields(category);
  }

  /**
   * Update subtype dropdown based on selected category
   */
  function updateSubtypeOptions(category) {
    storeOriginalSubtypes();

    // Remove all options except empty
    const options = Array.from(elements.subtypeField.querySelectorAll("option"));
    options.forEach((opt) => {
      if (opt.value) opt.remove();
    });

    // Add valid options for this category
    const validOptions = VALID_SUBTYPES[category] || [];
    validOptions.forEach((value) => {
      if (originalSubtypeOptions[value]) {
        const option = document.createElement("option");
        option.value = value;
        option.textContent = originalSubtypeOptions[value];
        elements.subtypeField.appendChild(option);
      }
    });

    // Auto-select if only one option
    if (validOptions.length === 1) {
      elements.subtypeField.value = validOptions[0];
    } else {
      elements.subtypeField.value = "";
    }
  }

  /**
   * Update owner field visibility and labels
   */
  function updateOwnerFields(category) {
    const categoryText = elements.categoryField.options[elements.categoryField.selectedIndex]?.text || "Owner";

    if (category === "p") {
      // Protein category
      elements.proteinOwnerField.style.display = "block";
      elements.ownerStateSelect.required = true;
      elements.otherOwnerField.style.display = "none";
      elements.ownerInput.required = false;
    } else {
      // Other categories
      elements.proteinOwnerField.style.display = "none";
      elements.ownerStateSelect.required = false;
      elements.otherOwnerField.style.display = "block";
      elements.ownerInput.required = true;
      elements.ownerTypeLabel.textContent = categoryText;
    }
  }

  /**
   * Update bio fields visibility (pH, solvent)
   */
  function updateBioFields(category) {
    const showBioFields = category === "d" || category === "p";

    if (showBioFields) {
      elements.bioFieldsRow.style.display = "";
    } else {
      elements.bioFieldsRow.style.display = "none";
      elements.phField.value = "";
      elements.solventField.value = "";
    }
  }

  /**
   * Update submit button text based on data presence
   */
  function updateSubmitButton() {
    const activeTab = document.querySelector("#data-source-tabs .nav-link.active");
    const isFileTab = activeTab?.id === "file-tab";
    const hasData = isFileTab ? elements.fileInput.value : elements.dataInput.value.trim();

    elements.submitBtn.textContent = hasData ? "Preview Spectrum" : "Submit";
  }

  /**
   * Handle form submission
   */
  function handleFormSubmit(e) {
    e.preventDefault();

    if (elements.previewSection.style.display !== "none") {
      // Preview already shown, submit the form
      submitFinalSpectrum();
    } else {
      // Show preview first
      showSpectrumPreview();
    }
  }

  /**
   * Show spectrum preview via AJAX
   */
  async function showSpectrumPreview() {
    const originalText = elements.submitBtn.textContent;
    elements.submitBtn.disabled = true;
    elements.submitBtn.textContent = "Processing...";

    try {
      // Prepare form data
      const formData = new FormData(elements.form);

      // Determine active tab and set data source
      const activeTab = document.querySelector("#data-source-tabs .nav-link.active");
      const isManualTab = activeTab?.id === "manual-tab";

      if (isManualTab) {
        formData.delete("file");
        formData.append("data_source", "manual");
      } else {
        formData.set("spectral_data", "");
        formData.append("data_source", "file");
      }

      // Get preview URL from form data attribute
      const previewUrl = elements.form.dataset.previewUrl;

      // Make request
      const response = await fetch(previewUrl, {
        method: "POST",
        headers: {
          "X-CSRFToken": getCsrfToken(),
          "X-Requested-With": "XMLHttpRequest",
        },
        body: formData,
      });

      const data = await response.json();

      if (response.ok && data.success) {
        currentPreviewData = data.preview;
        displaySpectrumPreview(data);
      } else {
        showError(data.error || "Failed to generate preview", data.details, data.form_errors);
      }
    } catch (error) {
      showError("Network error", error.message);
    } finally {
      elements.submitBtn.disabled = false;
      elements.submitBtn.textContent = originalText;
    }
  }

  /**
   * Display spectrum preview
   */
  function displaySpectrumPreview(response) {
    const preview = response.preview;

    // Clear any previous errors
    const alerts = elements.form.querySelectorAll(".alert-danger");
    alerts.forEach((alert) => alert.remove());

    // Update preview content
    elements.previewMessage.textContent = response.message;
    elements.previewPeakWave.textContent = preview.peak_wave || "N/A";
    elements.previewWaveRange.textContent = `${preview.min_wave}-${preview.max_wave}`;
    elements.previewDataPoints.textContent = preview.data_points;

    // Display SVG chart
    elements.previewChart.innerHTML = preview.svg;
    const svg = elements.previewChart.querySelector("svg");
    if (svg) {
      svg.style.maxWidth = "100%";
      svg.style.height = "auto";
      svg.style.display = "block";
      svg.style.margin = "0 auto";
    }

    // Show preview, hide form
    elements.formSection.style.display = "none";
    elements.previewSection.style.display = "block";

    // Scroll to preview
    elements.previewSection.scrollIntoView({ behavior: "smooth" });
  }

  /**
   * Hide preview and return to editing
   */
  function hidePreview() {
    if (!currentPreviewData) return;

    // Switch to manual tab
    const manualTab = document.getElementById("manual-tab");
    if (manualTab) {
      // Trigger Bootstrap tab
      if (window.$ && typeof window.$.fn.tab === "function") {
        window.$(manualTab).tab("show");
      } else {
        // Fallback without Bootstrap
        document.getElementById("file-tab")?.classList.remove("active");
        manualTab.classList.add("active");
        document.getElementById("file-panel")?.classList.remove("show", "active");
        document.getElementById("manual-panel")?.classList.add("show", "active");
      }
    }

    // Clear file input
    elements.fileInput.value = "";

    // Populate manual data field
    try {
      elements.dataInput.value = JSON.stringify(currentPreviewData.raw_data);
    } catch (e) {
      showError("Failed to load data for editing", e.message);
      return;
    }

    // Show form, hide preview
    elements.previewSection.style.display = "none";
    elements.formSection.style.display = "block";
    currentPreviewData = null;

    // Update button
    updateSubmitButton();

    // Focus on data field
    elements.dataInput.focus();
  }

  /**
   * Submit final spectrum
   */
  function submitFinalSpectrum() {
    if (!currentPreviewData) {
      showError("No preview data available");
      return;
    }

    // Add hidden field for data source
    const activeTab = document.querySelector("#data-source-tabs .nav-link.active");
    const dataSource = activeTab?.id === "manual-tab" ? "manual" : "file";

    const input = document.createElement("input");
    input.type = "hidden";
    input.name = "data_source";
    input.value = dataSource;
    elements.form.appendChild(input);

    // Submit form
    elements.form.submit();
  }

  /**
   * Show error message
   */
  function showError(message, details, formErrors) {
    let errorHtml = '<div class="alert alert-danger alert-dismissible fade show" role="alert">';
    errorHtml += `<strong>Error:</strong> ${message}`;

    if (details) {
      errorHtml += `<br><small>${details}</small>`;
    }

    if (formErrors && Object.keys(formErrors).length > 0) {
      errorHtml += "<br><br><strong>Form Issues:</strong><ul>";
      for (const [field, errors] of Object.entries(formErrors)) {
        errorHtml += `<li><strong>${field}:</strong> ${errors.join(", ")}</li>`;
      }
      errorHtml += "</ul>";
    }

    errorHtml += '<button type="button" class="close" data-dismiss="alert" aria-label="Close">';
    errorHtml += '<span aria-hidden="true">&times;</span></button></div>';

    // Insert at top of form
    elements.form.insertAdjacentHTML("afterbegin", errorHtml);

    // Scroll to error
    elements.form.scrollIntoView({ behavior: "smooth", block: "start" });
  }

  /**
   * Simple autocomplete for protein state field
   */
  function initStateAutocomplete() {
    const select = elements.ownerStateSelect;
    const searchUrl = elements.form.dataset.stateSearchUrl;

    if (!select || !searchUrl) return;

    // Add search input above select
    const wrapper = select.closest(".input-group") || select.parentElement;
    const searchInput = document.createElement("input");
    searchInput.type = "text";
    searchInput.className = "form-control mb-2";
    searchInput.placeholder = "Type to search proteins...";
    searchInput.id = "state-search-input";

    wrapper.insertBefore(searchInput, select);

    // Debounced search
    let searchTimeout;
    searchInput.addEventListener("input", (e) => {
      clearTimeout(searchTimeout);
      const query = e.target.value.trim();

      if (query.length < 2) {
        return;
      }

      searchTimeout = setTimeout(async () => {
        try {
          const response = await fetch(`${searchUrl}?q=${encodeURIComponent(query)}`, {
            headers: { "X-Requested-With": "XMLHttpRequest" },
          });
          const data = await response.json();

          // Update select options
          select.innerHTML = '<option value="">---------</option>';
          data.results.forEach((result) => {
            const option = document.createElement("option");
            option.value = result.id;
            option.textContent = result.text;
            select.appendChild(option);
          });
        } catch (error) {
          console.error("Search failed:", error);
        }
      }, 300);
    });
  }

  /**
   * Initialize all event listeners
   */
  function initEventListeners() {
    // Category change
    elements.categoryField.addEventListener("change", handleCategoryChange);

    // Form submission
    elements.form.addEventListener("submit", handleFormSubmit);

    // Data input changes
    elements.fileInput.addEventListener("change", updateSubmitButton);
    elements.dataInput.addEventListener("input", updateSubmitButton);

    // Tab changes
    document.querySelectorAll("#data-source-tabs a[data-toggle='tab']").forEach((tab) => {
      tab.addEventListener("shown.bs.tab", updateSubmitButton);
    });

    // Preview controls
    elements.editDataBtn.addEventListener("click", hidePreview);
    elements.submitFinalBtn.addEventListener("click", submitFinalSpectrum);

    // Clear state button
    if (elements.clearStateBtn) {
      elements.clearStateBtn.addEventListener("click", () => {
        elements.ownerStateSelect.value = "";
        const searchInput = document.getElementById("state-search-input");
        if (searchInput) searchInput.value = "";
      });
    }
  }

  /**
   * Initialize on DOM ready
   */
  function init() {
    initElements();
    initEventListeners();
    initStateAutocomplete();

    // Trigger initial category change to set up form state
    handleCategoryChange();
    updateSubmitButton();
  }

  // Run on DOM ready
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();

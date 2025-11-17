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
   * Modern autocomplete for protein state field
   * Replaces the select with a custom dropdown (similar to Select2 but lighter)
   */
  function initStateAutocomplete() {
    const select = elements.ownerStateSelect;
    const searchUrl = elements.form.dataset.stateSearchUrl;

    if (!select || !searchUrl) return;

    // Get the original input group or parent
    const originalWrapper = select.closest(".input-group");
    if (!originalWrapper) return;

    // Hide the original select but keep it for form submission
    select.style.display = "none";

    // Create custom autocomplete container
    const autocomplete = document.createElement("div");
    autocomplete.className = "custom-autocomplete";
    autocomplete.innerHTML = `
      <div class="autocomplete-input-wrapper">
        <input type="text"
               class="form-control autocomplete-input"
               placeholder="Type to search proteins..."
               id="state-search-input"
               autocomplete="off">
        <div class="autocomplete-spinner" style="display: none;">
          <span class="spinner-border spinner-border-sm text-primary"></span>
        </div>
      </div>
      <div class="autocomplete-dropdown" style="display: none;">
        <div class="autocomplete-results"></div>
      </div>
    `;

    // Insert before the input group append (clear button)
    const appendDiv = originalWrapper.querySelector(".input-group-append");
    originalWrapper.insertBefore(autocomplete, appendDiv);

    const searchInput = autocomplete.querySelector(".autocomplete-input");
    const dropdown = autocomplete.querySelector(".autocomplete-dropdown");
    const resultsContainer = autocomplete.querySelector(".autocomplete-results");
    const spinner = autocomplete.querySelector(".autocomplete-spinner");

    // Track current selection
    let selectedIndex = -1;
    let results = [];

    // Debounced search
    let searchTimeout;
    searchInput.addEventListener("input", async (e) => {
      clearTimeout(searchTimeout);
      const query = e.target.value.trim();

      if (query.length < 2) {
        hideDropdown();
        return;
      }

      spinner.style.display = "block";

      searchTimeout = setTimeout(async () => {
        try {
          const response = await fetch(`${searchUrl}?q=${encodeURIComponent(query)}`, {
            headers: { "X-Requested-With": "XMLHttpRequest" },
          });
          const data = await response.json();
          results = data.results || [];

          displayResults(results);
          spinner.style.display = "none";
        } catch (error) {
          console.error("Search failed:", error);
          spinner.style.display = "none";
          resultsContainer.innerHTML = '<div class="autocomplete-error">Search failed. Please try again.</div>';
          showDropdown();
        }
      }, 300);
    });

    // Display results
    function displayResults(items) {
      if (items.length === 0) {
        resultsContainer.innerHTML = '<div class="autocomplete-no-results">No proteins found</div>';
        showDropdown();
        return;
      }

      resultsContainer.innerHTML = items
        .map(
          (item, index) => `
        <div class="autocomplete-result" data-index="${index}" data-value="${item.id}">
          <div class="result-text">${item.text}</div>
        </div>
      `
        )
        .join("");

      // Add click handlers
      resultsContainer.querySelectorAll(".autocomplete-result").forEach((el) => {
        el.addEventListener("click", () => selectResult(parseInt(el.dataset.index)));
        el.addEventListener("mouseenter", () => {
          selectedIndex = parseInt(el.dataset.index);
          updateHighlight();
        });
      });

      selectedIndex = -1;
      showDropdown();
    }

    // Select a result
    function selectResult(index) {
      if (index < 0 || index >= results.length) return;

      const result = results[index];
      select.value = result.id;
      searchInput.value = result.text;
      hideDropdown();

      // Trigger change event on hidden select
      select.dispatchEvent(new Event("change", { bubbles: true }));
    }

    // Keyboard navigation
    searchInput.addEventListener("keydown", (e) => {
      if (!dropdown.style.display || dropdown.style.display === "none") return;

      switch (e.key) {
        case "ArrowDown":
          e.preventDefault();
          selectedIndex = Math.min(selectedIndex + 1, results.length - 1);
          updateHighlight();
          break;
        case "ArrowUp":
          e.preventDefault();
          selectedIndex = Math.max(selectedIndex - 1, -1);
          updateHighlight();
          break;
        case "Enter":
          e.preventDefault();
          if (selectedIndex >= 0) {
            selectResult(selectedIndex);
          }
          break;
        case "Escape":
          e.preventDefault();
          hideDropdown();
          break;
      }
    });

    // Update visual highlight
    function updateHighlight() {
      resultsContainer.querySelectorAll(".autocomplete-result").forEach((el, index) => {
        if (index === selectedIndex) {
          el.classList.add("active");
          el.scrollIntoView({ block: "nearest" });
        } else {
          el.classList.remove("active");
        }
      });
    }

    // Show/hide dropdown
    function showDropdown() {
      dropdown.style.display = "block";
    }

    function hideDropdown() {
      dropdown.style.display = "none";
      selectedIndex = -1;
    }

    // Close dropdown when clicking outside
    document.addEventListener("click", (e) => {
      if (!autocomplete.contains(e.target)) {
        hideDropdown();
      }
    });

    // Update clear button to also clear search input
    const clearBtn = elements.clearStateBtn;
    if (clearBtn) {
      clearBtn.addEventListener("click", () => {
        searchInput.value = "";
        hideDropdown();
      });
    }

    // If select already has a value, populate the search input
    if (select.value && select.selectedOptions[0]) {
      searchInput.value = select.selectedOptions[0].textContent;
    }
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

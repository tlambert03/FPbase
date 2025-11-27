import { fetchWithSentry } from "../ajax-sentry.js"

const VALIDATE_URL = "/ajax/validate_spectrumownername/"

/** Map subtype codes to display names (must match backend get_subtype_display) */
const SUBTYPE_DISPLAY_NAMES = {
  ex: "Excitation",
  ab: "Absorption",
  em: "Emission",
  "2p": "Two-Photon",
  bp: "Bandpass",
  bx: "Bandpass-Ex",
  bm: "Bandpass-Em",
  sp: "Shortpass",
  lp: "Longpass",
  bs: "Beamsplitter",
  qe: "Quantum Efficiency",
  pd: "Power Distribution",
}

/**
 * Check for similar spectrum owners and show warnings.
 *
 * @param {string} ownerName - The owner name to check
 * @param {string} category - The spectrum category (e.g., "D" for dye, "F" for filter)
 * @param {string} subtype - The spectrum subtype (e.g., "AB", "EM", "EX")
 * @param {HTMLElement} warningContainer - Element to display warnings
 * @returns {Promise<boolean>} True if exact match found (blocks submission), false otherwise
 */
export async function checkSimilarOwners(ownerName, category, subtype, warningContainer) {
  // Proteins use Select2 autocomplete - no duplicate checking needed
  if (!ownerName?.trim() || category === "P") {
    warningContainer.innerHTML = ""
    warningContainer.classList.add("d-none")
    return false
  }

  const formData = new URLSearchParams()
  formData.append("owner", ownerName)

  // Get CSRF token from form
  const csrfToken = document.querySelector("[name=csrfmiddlewaretoken]")?.value
  if (csrfToken) {
    formData.append("csrfmiddlewaretoken", csrfToken)
  }

  try {
    const response = await fetchWithSentry(VALIDATE_URL, {
      method: "POST",
      headers: {
        "Content-Type": "application/x-www-form-urlencoded",
        "X-Requested-With": "XMLHttpRequest",
      },
      body: formData,
    })

    const data = await response.json()

    if (!data.similars?.length) {
      warningContainer.innerHTML = ""
      warningContainer.classList.add("d-none")
      return false
    }

    // Filter similars based on category
    const relevantSimilars = filterSimilarsByCategory(data.similars, category, subtype)

    if (!relevantSimilars.length) {
      warningContainer.innerHTML = ""
      warningContainer.classList.add("d-none")
      return false
    }

    // Check for exact matches (case-insensitive) with same subtype
    // Backend returns display names ("Excitation"), we pass codes ("ex")
    const subtypeDisplayName = SUBTYPE_DISPLAY_NAMES[subtype.toLowerCase()]
    const exactMatches = relevantSimilars.filter((similar) => {
      const nameMatches = similar.name.toLowerCase() === ownerName.toLowerCase()
      if (!nameMatches) return false

      // Check if this exact owner already has this specific subtype
      return similar.spectra.includes(subtypeDisplayName)
    })

    if (exactMatches.length) {
      displayWarning(exactMatches, warningContainer, true)
      return true
    }

    displayWarning(relevantSimilars, warningContainer, false)
    return false
  } catch (error) {
    console.error("Error checking similar owners:", error)
    return false
  }
}

/**
 * Filter similar owners based on category and subtype.
 *
 * For dyes: only show if the similar dye has the SAME subtype
 * For others: show all similar names (they only have one spectrum per owner)
 */
function filterSimilarsByCategory(similars, category, subtype) {
  if (category === "D") {
    // Dye: only warn if similar dye has the same subtype
    // Backend returns display names, so convert code to display name
    const subtypeDisplayName = SUBTYPE_DISPLAY_NAMES[subtype.toLowerCase()]
    return similars.filter((similar) => similar.spectra.includes(subtypeDisplayName))
  }
  // Filter/Camera/Light: warn for any similar names
  return similars
}

/**
 * Display warning message with links to similar owners.
 *
 * @param {Array} similars - Array of similar spectrum owners
 * @param {HTMLElement} warningContainer - Element to display warnings
 * @param {boolean} isExactMatch - True if this is an exact match (error), false for similar matches (warning)
 */
function displayWarning(similars, warningContainer, isExactMatch) {
  const items = similars
    .map((similar) => {
      const spectraText = similar.spectra.length ? ` (${similar.spectra.join(", ")})` : ""
      const spectrumIds = similar.spectrum_ids?.join(",") || ""
      const url = spectrumIds ? `/spectra/?s=${spectrumIds}` : similar.url
      return `<li><a href="${url}" target="_blank" class="text-danger">${similar.name}</a>${spectraText}</li>`
    })
    .join("")

  if (isExactMatch) {
    warningContainer.className = "alert alert-danger small mt-2 mb-0"
    warningContainer.innerHTML = `
      <strong>🚫 Exact match found:</strong> This spectrum already exists:
      <ul class="mb-0 mt-1">${items}</ul>
      <div class="small mt-1"><em>Cannot submit duplicate spectrum</em></div>
    `
  } else {
    warningContainer.className = "alert alert-warning small mt-2 mb-0"
    warningContainer.innerHTML = `
      <strong>⚠️ Possible duplicate:</strong> Similarly named existing spectra:
      <ul class="mb-0 mt-1">${items}</ul>
    `
  }

  warningContainer.classList.remove("d-none")
}

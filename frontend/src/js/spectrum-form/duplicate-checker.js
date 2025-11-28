import { fetchWithSentry } from "../ajax-sentry.js"

const VALIDATE_URL = "/ajax/validate_spectrumownername/"

/** Map subtype codes to display names (must match backend Spectrum.SUBTYPE_CHOICES) */
const SUBTYPE_DISPLAY_NAMES = {
  ex: "Excitation",
  ab: "Absorption",
  em: "Emission",
  "2p": "Two Photon Abs",
  bp: "Bandpass",
  bx: "Bandpass-Ex",
  bm: "Bandpass-Em",
  sp: "Shortpass",
  lp: "Longpass",
  bs: "Beamsplitter",
  qe: "Quantum Efficiency",
  pd: "Power Distribution",
}

/** Map display names back to subtype codes */
const DISPLAY_NAME_TO_CODE = Object.fromEntries(
  Object.entries(SUBTYPE_DISPLAY_NAMES).map(([code, name]) => [name, code])
)

/**
 * Check for similar spectrum owners, show warnings, and return existing subtypes.
 *
 * @param {string} ownerName - The owner name to check
 * @param {string} category - The spectrum category code (d=dye, p=protein, f=filter, c=camera, l=light)
 * @param {HTMLElement} warningContainer - Element to display warnings
 * @param {Object} [options] - Additional options
 * @param {boolean} [options.blockOnExactMatch] - If true, show error (red) for exact matches
 * @param {boolean} [options.suppressWarningOnExactMatch] - If true, hide warning when exact match found
 *   (subtypes still returned for disabling options)
 * @returns {Promise<{hasExactMatch: boolean, existingSubtypes: string[]}>}
 *   - hasExactMatch: True if exact name match found (for any subtype)
 *   - existingSubtypes: Array of subtype codes that already exist for the exact match
 */
export async function checkSimilarOwners(ownerName, category, warningContainer, options = {}) {
  const { blockOnExactMatch = false, suppressWarningOnExactMatch = false } = options
  const result = { hasExactMatch: false, existingSubtypes: [] }

  if (!ownerName?.trim()) {
    warningContainer.innerHTML = ""
    warningContainer.classList.add("d-none")
    return result
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
      return result
    }

    // Filter similars based on category
    const relevantSimilars = filterSimilarsByCategory(data.similars, category)

    if (!relevantSimilars.length) {
      warningContainer.innerHTML = ""
      warningContainer.classList.add("d-none")
      return result
    }

    // Find exact name matches (case-insensitive)
    const exactMatch = relevantSimilars.find(
      (similar) => similar.name.toLowerCase() === ownerName.toLowerCase()
    )

    if (exactMatch) {
      result.hasExactMatch = true
      // Convert display names back to codes
      result.existingSubtypes = exactMatch.spectra
        .map((displayName) => DISPLAY_NAME_TO_CODE[displayName])
        .filter(Boolean)
    }

    // Suppress warning for exact matches if requested (e.g., proteins where we just disable subtypes)
    if (suppressWarningOnExactMatch && result.hasExactMatch) {
      // Only show exact match, no "similar" warning needed
      warningContainer.innerHTML = ""
      warningContainer.classList.add("d-none")
      return result
    }

    // Show warning or error based on options
    const showAsError = blockOnExactMatch && result.hasExactMatch
    displayWarning(relevantSimilars, warningContainer, showAsError)
    return result
  } catch (error) {
    console.error("Error checking similar owners:", error)
    return result
  }
}

/**
 * Filter similar owners to only show those matching the selected category.
 *
 * This ensures users only see potential duplicates within the same category
 * (e.g., when entering a dye, only show similar dyes, not proteins or filters).
 */
function filterSimilarsByCategory(similars, category) {
  return similars.filter((similar) => similar.category === category)
}

/**
 * Display warning message with links to similar owners.
 *
 * @param {Array} similars - Array of similar spectrum owners
 * @param {HTMLElement} warningContainer - Element to display warnings
 * @param {boolean} [isBlockingError] - If true, show as error (red) instead of warning (yellow)
 */
function displayWarning(similars, warningContainer, isBlockingError = false) {
  const items = similars
    .map((similar) => {
      const spectraText = similar.spectra.length ? ` (${similar.spectra.join(", ")})` : ""
      const spectrumIds = similar.spectrum_ids?.join(",") || ""
      const url = spectrumIds ? `/spectra/?s=${spectrumIds}` : similar.url
      return `<li><a href="${url}" target="_blank" class="alert-link">${similar.name}</a>${spectraText}</li>`
    })
    .join("")

  if (isBlockingError) {
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

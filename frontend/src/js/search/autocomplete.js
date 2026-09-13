import { escapeHtml, highlight, SearchIndex } from "./engine.js"

const $ = window.jQuery // jQuery loaded from CDN

const PROTEIN_FIELDS = [
  { key: "name", weight: 1 },
  { key: "aliases", weight: 0.85 },
  { key: "uuid", weight: 0.9, id: true },
  { key: "pdb", weight: 0.9, id: true },
  { key: "genbank", weight: 0.9, id: true },
  { key: "uniprot", weight: 0.9, id: true },
  { key: "ipg_id", weight: 0.9, id: true },
]
const REFERENCE_FIELDS = [
  { key: "primary", weight: 1 },
  { key: "citation", weight: 0.95 },
  { key: "title", weight: 0.7 },
  { key: "doi", weight: 0.9, id: true },
  { key: "pmid", weight: 0.9, id: true },
  { key: "secondary", weight: 0.6 },
]
// organisms only match whole words by prefix: substrings and short typos are noise here
const ORGANISM_OPTIONS = { minInfix: Infinity, minTypo1: 5, fallback: false }

const PROTEIN_HINT_LABELS = {
  aliases: "aka",
  uuid: "FPbase ID",
  pdb: "PDB",
  genbank: "GenBank",
  uniprot: "UniProt",
  ipg_id: "IPG ID",
}
const REFERENCE_HINT_LABELS = {
  doi: "DOI",
  pmid: "PMID",
  primary: "Protein",
  secondary: "2˚ Protein",
}
const SWITCH_ABBREVIATIONS = {
  photoswitchable: "PS",
  photoactivatable: "PA",
  photoconvertible: "PC",
  "multi-photochromic": "MPC",
  multistate: "MS",
  timer: "Time",
}

let indexesPromise = null
let indexes = null

/** Fetch the search index (browser-cached, see /api/search-index/) and build the indexes. */
function loadIndexes() {
  indexesPromise ||= fetch(window.FPBASE.searchIndexURL)
    .then((response) => {
      if (!response.ok) throw new Error(`Search index request failed: ${response.status}`)
      return response.json()
    })
    .then((data) => {
      indexes = {
        proteins: new SearchIndex(data.proteins, PROTEIN_FIELDS),
        references: new SearchIndex(data.references, REFERENCE_FIELDS),
        organisms: new SearchIndex(data.organisms, [{ key: "name", weight: 1 }], ORGANISM_OPTIONS),
      }
    })
    .catch((error) => {
      indexesPromise = null // allow a retry on the next keystroke
      window.Sentry?.captureException(error)
      throw error
    })
  return indexesPromise
}

// autocomplete.js JSON-serializes every suggestion into the DOM, so suggestions are
// rendered here and only small {url, display, html} objects are handed over.
function source(name, limit, render, displayKey) {
  const run = (query, callback) => {
    const index = indexes[name]
    callback(
      index.search(query, limit).map((hit) => ({
        url: hit.record.url,
        display: hit.record[displayKey],
        html: render({ ...hit, index, query }),
      }))
    )
  }
  return (query, callback) => {
    // answer synchronously once loaded: if this dataset rendered after the (synchronous)
    // footer, autoselect would put the cursor on "Advanced search" and Enter would go there
    if (indexes) return run(query, callback)
    loadIndexes()
      .then(() => run(query, callback))
      .catch(() => callback([]))
  }
}

function hints(hit, exclude = []) {
  const byKey = {}
  for (const { key, value } of hit.index.matchedValues(hit, hit.query, exclude)) {
    byKey[key] ??= []
    byKey[key].push(highlight(value, hit.query))
  }
  return byKey
}

function proteinSuggestion(hit) {
  const p = hit.record
  let col = "gray50"
  if (p.switch && p.switch !== "Basic") {
    col = "rainbow"
  } else if (p.color && !p.color.includes("Stokes")) {
    col = p.color.toLowerCase().replace(/ |\//g, "_")
  }
  let str = `<img class='type protein' src='${window.FPBASE.imageDir}gfp_${col}_40.png'>`
  str += highlight(p.name, hit.query)
  if (p.spectra) {
    const src = `/spectra_img/${encodeURIComponent(p.slug)}.png?xlabels=0&xlim=400,800`
    str += `<img class='spectra' src='${src}'>`
  }
  const matched = hints(hit, ["name"])
  const parts = Object.entries(matched).map(
    ([key, values]) => `${PROTEIN_HINT_LABELS[key]}: ${values.join(", ")}`
  )
  if (parts.length) str += `<span class='highlighted-hits'>(${parts.join("; ")})</span>`
  let info = ""
  if (p.switch === "Basic") {
    if (p.ex && p.em) info = `${p.ex}/${p.em}`
  } else if (p.switch) {
    info = SWITCH_ABBREVIATIONS[p.switch.toLowerCase()] || ""
  }
  str += `<span class='info'>${escapeHtml(info)}</span>`
  return `<a href='${escapeHtml(p.url)}'><div>${str}</div></a>`
}

function referenceSuggestion(hit) {
  const ref = hit.record
  let str = highlight(ref.citation, hit.query)
  str += `<img class='type' src='${window.FPBASE.imageDir}ref.png'>`
  const matched = hints(hit, ["citation"])
  if (matched.primary) delete matched.secondary
  const parts = Object.keys(REFERENCE_HINT_LABELS)
    .filter((key) => matched[key])
    .map((key) => `${REFERENCE_HINT_LABELS[key]}: ${matched[key].join(", ")}`)
  if (parts.length) str += `<span class='highlighted-hits'>(${parts.join("; ")})</span>`
  if (matched.title) str += `<div class="ref-title">${matched.title[0]}</div>`
  return `<a href='${escapeHtml(ref.url)}'><div>${str}</div></a>`
}

function organismSuggestion(hit) {
  const org = hit.record
  const str = `${highlight(org.name, hit.query)}<img class='type' src='${window.FPBASE.imageDir}organism_icon.png'>`
  return `<a href='${escapeHtml(org.url)}'><div>${str}</div></a>`
}

function empty({ query }) {
  const href = query ? `/search/?name__icontains=${encodeURIComponent(query.trim())}` : "/search/"
  return `<div class="empty"><span class="nohits"></span>No results... try the <a href="${escapeHtml(href)}">advanced search</a></div>`
}

// Guard to prevent double initialization
let isInitialized = false

/**
 * Wait for autocomplete.js library to be available
 * @param {number} maxWaitMs - Maximum time to wait in milliseconds
 * @returns {Promise<boolean>} - Resolves to true if available, false if timeout
 */
async function waitForAutocomplete(maxWaitMs = 2000) {
  const startTime = Date.now()
  const checkInterval = 50

  while (Date.now() - startTime < maxWaitMs) {
    if (typeof $.fn.autocomplete !== "undefined") {
      return true
    }
    await new Promise((resolve) => setTimeout(resolve, checkInterval))
  }
  return false
}

/**
 * Initialize the site search autocomplete
 * Must be called after DOM is ready and autocomplete.js is loaded
 */
export default async function initAutocomplete() {
  // Prevent double initialization
  if (isInitialized) {
    return
  }

  // Wait for search input to be available in DOM
  const $searchInput = $("#algolia-search-input")
  if (!$searchInput.length) {
    console.warn("Search input not found in DOM")
    return
  }

  // start downloading the index as soon as the user shows interest in searching
  $searchInput.one("focus pointerenter touchstart", () => loadIndexes().catch(() => {}))
  if (document.activeElement === $searchInput[0]) loadIndexes().catch(() => {})

  // Wait for autocomplete.js library (loaded from CDN with defer)
  const autocompleteAvailable = await waitForAutocomplete()
  if (!autocompleteAvailable) {
    console.error("Autocomplete plugin failed to load after 2 seconds")
    if (window.Sentry) {
      Sentry.captureMessage("Autocomplete CDN script failed to load", "warning")
    }
    return
  }

  isInitialized = true

  // Initialize autocomplete on the search input
  $searchInput
    .autocomplete(
      {
        minLength: 3,
        autoselect: true,
        autoselectOnBlur: window.mobilecheck(),
        templates: {
          empty: empty,
        },
      },
      [
        {
          source: source("proteins", 5, proteinSuggestion, "name"),
          displayKey: "display",
          templates: { suggestion: (suggestion) => suggestion.html },
        },
        {
          source: source("references", 3, referenceSuggestion, "citation"),
          displayKey: "display",
          templates: { suggestion: (suggestion) => suggestion.html },
        },
        {
          source: source("organisms", 2, organismSuggestion, "name"),
          displayKey: "display",
          templates: { suggestion: (suggestion) => suggestion.html },
        },
        {
          source: (query, callback) => {
            const footer = () =>
              callback([{ query, url: `/search/?q=${encodeURIComponent(query)}` }])
            // render after the other datasets so autoselect picks a real result
            if (indexes) footer()
            else loadIndexes().then(footer, footer)
          },
          templates: {
            suggestion: (suggestion) =>
              `<div class="search-footer"><a class="asearch" href="${escapeHtml(suggestion.url)}">` +
              `Advanced search for: <em>${escapeHtml(suggestion.query)}</em></a></div>`,
          },
        },
      ]
    )
    .on("autocomplete:selected", (_event, suggestion, _dataset, context) => {
      if (context.selectionMethod === "click") {
        return
      }
      // Change the page, for example, on other events
      window.location.assign(suggestion.url)
    })

  const $hintInput = $searchInput.parent().find(".aa-hint")
  if ($hintInput.length) {
    $hintInput.attr("name", "search-hint")
    $hintInput.attr("id", "algolia-search-hint")
  }
}

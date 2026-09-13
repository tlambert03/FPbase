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
  { key: "title", weight: 0.7, typos: false },
  { key: "doi", weight: 0.9, id: true },
  { key: "pmid", weight: 0.9, id: true },
  { key: "secondary", weight: 0.6 },
]
// organisms only match whole words by prefix: substrings and short typos are noise here
const ORGANISM_OPTIONS = { minInfix: Infinity, minTypo1: 5, fallback: false }
// abbreviations whose expansion isn't already part of the dye names (JF646, BUV395 are)
const DYE_OPTIONS = {
  synonyms: { af: "alexa fluor", bv: "brilliant violet", sb: "super bright" },
}

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

// Search analytics: one GA4 `search` event per search, sent when a result is picked, the
// advanced search is used, or the search is abandoned (blur, or leaving the page).
let typedQuery = "" // the query as typed (arrow keys replace the input's value)
let tracked = false
const resultCounts = {} // dataset -> {query, n}

function trackSearch(params) {
  if (tracked || typedQuery.length < 3) return
  tracked = true
  const count = Object.values(resultCounts)
    .filter((c) => c.query === typedQuery)
    .reduce((n, c) => n + c.n, 0)
  window.gtag?.("event", "search", {
    search_term: typedQuery,
    result_count: count,
    transport_type: "beacon",
    ...params,
  })
}

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
        dyes: new SearchIndex(data.dyes, [{ key: "name", weight: 1 }], DYE_OPTIONS),
      }
    })
    .catch((error) => {
      indexesPromise = null // allow a retry on the next keystroke
      window.Sentry?.captureException(error)
      throw error
    })
  return indexesPromise
}

const LIMITS = { proteins: 5, dyes: 3, references: 3, organisms: 2 }
let lastResults = { query: null, results: {} }

/** Search every section at once, hiding sections that only have typo matches
 * when another section has a real match ("cy5": the dye, not CyPet). */
function resultsFor(query) {
  if (lastResults.query !== query) {
    const results = {}
    for (const [name, limit] of Object.entries(LIMITS)) {
      results[name] = indexes[name].search(query, limit)
    }
    const hits = Object.values(results)
    if (hits.some((section) => section.some((h) => !h.fuzzy))) {
      for (const [name, section] of Object.entries(results)) {
        if (section.every((h) => h.fuzzy)) results[name] = []
      }
    }
    lastResults = { query, results }
  }
  return lastResults.results
}

// autocomplete.js JSON-serializes every suggestion into the DOM, so suggestions are
// rendered here and only small objects are handed over.
function source(name, type, render, displayKey) {
  const run = (query, callback) => {
    const index = indexes[name]
    const hits = resultsFor(query)[name]
    typedQuery = query.trim()
    resultCounts[name] = { query: typedQuery, n: hits.length }
    callback(
      hits.map((hit, i) => ({
        url: hit.record.url,
        display: hit.record[displayKey],
        html: render({ ...hit, index }),
        type,
        rank: i + 1,
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
  for (const { key, value } of hit.index.matchedValues(hit, exclude)) {
    byKey[key] ??= []
    byKey[key].push(highlight(value, hit.tokens))
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
  str += highlight(p.name, hit.tokens)
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
  let str = highlight(ref.citation, hit.tokens)
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
  const str = `${highlight(org.name, hit.tokens)}<img class='type' src='${window.FPBASE.imageDir}organism_icon.png'>`
  return `<a href='${escapeHtml(org.url)}'><div>${str}</div></a>`
}

// a small aromatic ring, filled with the dye's emission color
function dyeIcon(color) {
  const fill = /^#[0-9a-f]{3,6}$/i.test(color ?? "") ? color : "#bbb"
  return (
    `<svg class='type dye' viewBox='0 0 24 24' aria-hidden='true'>` +
    `<path d='M12 2.5 20.2 7.25v9.5L12 21.5 3.8 16.75v-9.5Z' fill='${fill}' stroke='#444' stroke-opacity='.6' stroke-width='1.5'/>` +
    `<circle cx='12' cy='12' r='4.2' fill='none' stroke='#444' stroke-opacity='.6' stroke-width='1.3'/></svg>`
  )
}

function dyeSuggestion(hit) {
  const dye = hit.record
  let str = dyeIcon(dye.color) + highlight(dye.name, hit.tokens)
  if (dye.ex && dye.em) str += `<span class='info'>${escapeHtml(`${dye.ex}/${dye.em}`)}</span>`
  return `<a href='${escapeHtml(dye.url)}'><div>${str}</div></a>`
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
        // sections without hits render nothing, so e.g. a dye search shows only dyes
        {
          source: source("proteins", "protein", proteinSuggestion, "name"),
          displayKey: "display",
          templates: { suggestion: (suggestion) => suggestion.html },
        },
        {
          source: source("dyes", "dye", dyeSuggestion, "name"),
          displayKey: "display",
          templates: { suggestion: (suggestion) => suggestion.html },
        },
        {
          source: source("references", "reference", referenceSuggestion, "citation"),
          displayKey: "display",
          templates: { suggestion: (suggestion) => suggestion.html },
        },
        {
          source: source("organisms", "organism", organismSuggestion, "name"),
          displayKey: "display",
          templates: { suggestion: (suggestion) => suggestion.html },
        },
        {
          source: (query, callback) => {
            const footer = () =>
              callback([
                { query, url: `/search/?q=${encodeURIComponent(query)}`, type: "advanced" },
              ])
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
      trackSearch({ result_type: suggestion.type, result_rank: suggestion.rank ?? 0 })
      if (context.selectionMethod === "click") {
        return
      }
      // Change the page, for example, on other events
      window.location.assign(suggestion.url)
    })

  $searchInput.on("input", () => {
    tracked = false
    if ($searchInput.val().trim().length < 3) typedQuery = ""
  })
  $searchInput.on("blur", () => trackSearch({ result_type: "none" }))
  // Enter without a highlighted suggestion submits the form to the advanced search
  $searchInput.closest("form").on("submit", () => trackSearch({ result_type: "advanced" }))
  window.addEventListener("pagehide", () => trackSearch({ result_type: "none" }))

  const $hintInput = $searchInput.parent().find(".aa-hint")
  if ($hintInput.length) {
    $hintInput.attr("name", "search-hint")
    $hintInput.attr("id", "algolia-search-hint")
  }
}

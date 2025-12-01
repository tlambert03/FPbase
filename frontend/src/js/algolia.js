import { autocomplete } from "@algolia/autocomplete-js"
import "@algolia/autocomplete-theme-classic"

// Algolia has a 512 byte limit for query strings
// Use 500 to be safe (accounting for multi-byte characters)
const MAX_QUERY_LENGTH = 500

function checkObject(val, prop, str) {
  const propDict = {
    genbank: "GenBank",
    pdb: "PDB",
    uuid: "FPbase ID",
    ipg_id: "IPG ID",
    uniprot: "UniProt",
    aliases: "aka",
    doi: "DOI",
    pmid: "PMID",
    title: "Title",
    _excerpts: "Excerpt",
    prot_primary: "Protein",
    prot_secondary: "Protein",
  }
  if (val.matchLevel === "full" || val.matchLevel === "partial") {
    if (str.length > 0) {
      str = `${str}; `
    }
    str = `${str + propDict[prop]}: ${val.value}`
  }
  return str
}

function highlightHits(high) {
  let str = ""
  for (const prop in high) {
    if (Object.hasOwn(high, prop) && prop !== "name" && prop !== "citation") {
      if (high[prop].constructor === Array) {
        for (let i = 0; i < high[prop].length; i++) {
          if (typeof high[prop][i] === "object") {
            str = checkObject(high[prop][i], prop, str)
          }
        }
      } else {
        if (typeof high[prop] === "object") {
          str = checkObject(high[prop], prop, str)
        }
      }
    }
  }
  if (str) {
    return `<span class='highlighted-hits'>(${str})</span>`
  } else {
    return ""
  }
}

function highlightRefHits(high) {
  function recurseMatches(obj) {
    const results = {}

    function innerRecurse(obj, _key) {
      if (_key !== undefined) {
        obj = obj[_key]
      }
      for (const key in obj) {
        if (Object.hasOwn(obj, key)) {
          if (obj[key].constructor === Array) {
            innerRecurse(obj, key)
          } else if (Object.hasOwn(obj[key], "matchLevel")) {
            if (obj[key].matchLevel !== "none" && obj[key].value !== undefined) {
              let attr = key
              if (_key !== undefined) {
                attr = _key
              }
              if (!Object.hasOwn(results, attr)) {
                results[attr] = []
              }
              results[attr].push([obj[key].value, obj[key].matchLevel])
            }
          }
        }
      }
    }
    innerRecurse(obj)
    return results
  }

  function truncate(str, no_words) {
    return `${str.split(" ").splice(0, no_words).join(" ")} ...`
  }

  const results = recurseMatches(high)

  if (Object.hasOwn(results, "prot_primary")) {
    delete results.prot_secondary
  }
  let str = ""
  const items = [
    ["doi", "DOI"],
    ["pmid", "PMID"],
    ["prot_primary", "Protein"],
    ["prot_secondary", "2˚ Protein"],
  ]
  for (let x = 0; x < items.length; x++) {
    const key = items[x][0]
    const title = items[x][1]
    if (Object.hasOwn(results, key)) {
      const _str = []
      for (let i = 0; i < results[key].length; i++) {
        _str.push(results[key][i][0])
      }
      str = `${str + title}: ${_str.join(", ")}`
    }
  }
  if (str) {
    str = `<span class='highlighted-hits'>(${str})</span>`
  }
  if (Object.hasOwn(results, "title")) {
    str = `${str}<div class="ref-title" >${results.title[0][0]}</div>`
  } else if (Object.hasOwn(results, "_excerpts")) {
    if (results._excerpts.some((d) => d[1] === "full")) {
      for (let e = 0; e < results._excerpts.length; e++) {
        if (results._excerpts[e][1] === "full") {
          const exc = results._excerpts[e][0]
          const pre = exc.split("<em>")[0].split(" ")
          let _pre = ""
          if (pre.length > 5) {
            _pre = `${_pre}... `
          }
          str =
            str +
            '<div class="excerpt" >"' +
            _pre +
            pre.slice(Math.max(pre.length - 5, 0)).join(" ") +
            "<em>" +
            truncate(exc.split("<em>")[1], 7) +
            '"</div>'
          break
        }
      }
    }
  }
  return str
}

// Guard to prevent double initialization
let isInitialized = false

/**
 * Initialize Algolia autocomplete search
 * Must be called after DOM is ready
 */
export default async function initAutocomplete() {
  // Prevent double initialization
  if (isInitialized) {
    return
  }

  // Wait for search input to be available in DOM
  const searchInput = document.getElementById("algolia-search-input")
  if (!searchInput) {
    console.warn("Algolia search input not found in DOM")
    return
  }

  // Mark as initialized before async import
  isInitialized = true

  const algoliasearch = await import("algoliasearch/lite")

  const searchClient = algoliasearch.default(
    window.FPBASE.ALGOLIA.appID,
    window.FPBASE.ALGOLIA.publicKey
  )

  // Create a container div and replace the input
  const container = document.createElement("div")
  container.id = "autocomplete-container"
  searchInput.parentNode.replaceChild(container, searchInput)

  // Initialize autocomplete
  autocomplete({
    container: "#autocomplete-container",
    placeholder: "Search",
    autoselect: true,
    openOnFocus: false,
    detachedMediaQuery: "none", // Keep dropdown attached
    getSources({ query }) {
      // Truncate query if it exceeds Algolia's limit
      const truncatedQuery =
        query.length > MAX_QUERY_LENGTH ? query.substring(0, MAX_QUERY_LENGTH) : query

      if (truncatedQuery.length < 3) {
        return []
      }

      return [
        // Proteins source
        {
          sourceId: "proteins",
          async getItems() {
            const { results } = await searchClient.search([
              {
                indexName: window.FPBASE.ALGOLIA.proteinIndex,
                query: truncatedQuery,
                hitsPerPage: 5,
              },
            ])
            return results[0].hits
          },
          templates: {
            item({ item, html }) {
              let col
              if (item.switchType && item.switchType !== "Basic") {
                col = "rainbow"
              } else if (item.color && !item.color.includes("Stokes")) {
                col = item.color.toLowerCase().replace(/ |\//g, "_")
              } else {
                col = "gray50"
              }
              const imgSrc = `${window.FPBASE.imageDir}gfp_${col}_40.png`
              let content = `<img class='type protein' src='${imgSrc}'>`
              content += item._highlightResult.name.value
              if (item.img_url) {
                content += `<img class='spectra' src='${item.img_url}'>`
              }
              content += highlightHits(item._highlightResult)
              let info = ""
              if (item.switchType === "Basic") {
                if (item.ex && item.em) {
                  info = `${item.ex}/${item.em}`
                }
              } else if (item.switchType) {
                info =
                  {
                    photoswitchable: "PS",
                    photoactivatable: "PA",
                    photoconvertible: "PC",
                    "multi-photochromic": "MPC",
                    multistate: "MS",
                    timer: "Time",
                  }[item.switchType.toLowerCase()] || ""
              }
              content += `<span class='info'>${info}</span>`
              return html`<a href="${item.url}"><div>${html([content])}</div></a>`
            },
          },
          onSelect({ item }) {
            window.location.assign(item.url)
          },
        },
        // References source
        {
          sourceId: "references",
          async getItems() {
            const { results } = await searchClient.search([
              {
                indexName: window.FPBASE.ALGOLIA.referenceIndex,
                query: truncatedQuery,
                hitsPerPage: 3,
              },
            ])
            return results[0].hits
          },
          templates: {
            item({ item, html }) {
              let content = item._highlightResult.citation.value
              content += `<img class='type' src='${window.FPBASE.imageDir}ref.png'>`
              content += highlightRefHits(item._highlightResult)
              return html`<a href="${item.url}"><div>${html([content])}</div></a>`
            },
          },
          onSelect({ item }) {
            window.location.assign(item.url)
          },
        },
        // Organisms source
        {
          sourceId: "organisms",
          async getItems() {
            const { results } = await searchClient.search([
              {
                indexName: window.FPBASE.ALGOLIA.organismIndex,
                query: truncatedQuery,
                hitsPerPage: 2,
              },
            ])
            return results[0].hits
          },
          templates: {
            item({ item, html }) {
              let content = item._highlightResult.scientific_name.value
              content += `<img class='type' src='${window.FPBASE.imageDir}organism_icon.png'>`
              return html`<a href="${item.url}"><div>${html([content])}</div></a>`
            },
          },
          onSelect({ item }) {
            window.location.assign(item.url)
          },
        },
        // Advanced search footer
        {
          sourceId: "advanced-search",
          getItems() {
            return [
              {
                query: truncatedQuery,
                url: `/search/?q=${encodeURI(truncatedQuery)}`,
              },
            ]
          },
          templates: {
            item({ item, html }) {
              const content =
                '<div class="search-footer"><a class="asearch" href="' +
                item.url +
                '">Advanced search for: <em>' +
                item.query +
                '</em></a><div class="branding">search powered by <a href="https://algolia.com"><img src="' +
                window.FPBASE.imageDir +
                'logo-algolia-nebula-blue-full.svg" /></a></div></div>'
              return html([content])
            },
          },
          onSelect({ item }) {
            window.location.assign(item.url)
          },
        },
      ]
    },
  })
}

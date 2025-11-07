import { autocomplete } from "@algolia/autocomplete-js"
import { createLocalStorageRecentSearchesPlugin } from "@algolia/autocomplete-plugin-recent-searches"
import { liteClient } from "algoliasearch/lite"
import { useEffect, useRef } from "react"
import "@algolia/autocomplete-theme-classic"
import "./SearchAutocomplete.css"

// ============================================================================
// Constants
// ============================================================================

const IMAGE_DIR = () => window.FPBASE?.imageDir || "/static/images/"

const SEARCH_CONFIG = {
  placeholder: "Search proteins, organisms, references...",
  recentSearchesKey: "fpbase-search",
  recentSearchesLimit: 5,
  detachedMediaQuery: "(max-width: 768px)",
  debugQuery: "GFP",
}

const SOURCE_CONFIG = {
  proteins: {
    hitsPerPage: 5,
    attributes: ["name", "uuid", "aliases", "switchType", "color", "ex", "em", "img_url", "url"],
    highlightAttributes: ["name", "aliases"],
  },
  references: {
    hitsPerPage: 3,
    highlightAttributes: ["citation", "title"],
  },
  organisms: {
    hitsPerPage: 2,
    highlightAttributes: ["scientific_name"],
  },
}

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Determines the icon color for a protein based on its properties
 * @param {Object} protein - The protein hit object
 * @returns {string} The icon color identifier
 */
function getProteinIconColor(protein) {
  if (protein.switchType && protein.switchType !== "Basic") {
    return "rainbow"
  }
  if (protein.color && !protein.color.includes("Stokes")) {
    return protein.color.toLowerCase().replace(/[ /]/g, "_")
  }
  return "gray50"
}

/**
 * Filters aliases to only those that contain search matches (highlighted with <em> tags)
 * @param {Object} hit - The Algolia hit object
 * @returns {Array} Array of matched alias highlight results
 */
function getMatchedAliases(hit) {
  if (!hit._highlightResult?.aliases) return []
  return hit._highlightResult.aliases.filter((alias) => alias.value.includes("<em>"))
}

/**
 * Gets the highlighted or fallback value for a field
 * @param {Object} hit - The Algolia hit object
 * @param {string} field - The field name
 * @returns {string} The highlighted or original value
 */
function getHighlightedValue(hit, field) {
  return hit._highlightResult?.[field]?.value || hit[field]
}

/**
 * Formats aliases for display
 * @param {Array} aliases - Array of alias highlight results
 * @returns {string} Formatted alias string
 */
function formatAliases(aliases) {
  return (
    "aka: " +
    aliases
      .slice(0, 3)
      .map((a) => a.value)
      .join(", ")
  )
}

// ============================================================================
// Hit Template Components
// ============================================================================

/**
 * Renders a protein search result
 * @param {Object} params
 * @param {Object} params.hit - The protein hit data
 * @param {Function} params.html - Algolia's HTML template function
 */
function ProteinHit({ hit, html }) {
  const imageDir = IMAGE_DIR()
  const iconColor = getProteinIconColor(hit)
  const matchedAliases = getMatchedAliases(hit)

  return html`
    <a href="${hit.url}" class="aa-ItemLink">
      <div class="aa-ItemContent aa-ItemContent--protein">
        <img
          src="${imageDir}gfp_${iconColor}_40.png"
          alt="Protein icon"
          class="aa-ItemIcon aa-ItemIcon--protein"
          loading="lazy"
        />
        <div class="aa-ItemBody">
          <div
            class="aa-ItemTitle"
            dangerouslySetInnerHTML=${{ __html: getHighlightedValue(hit, "name") }}
          ></div>
          ${
            matchedAliases.length > 0
              ? html`<div
                class="aa-ItemAliases"
                dangerouslySetInnerHTML=${{ __html: formatAliases(matchedAliases) }}
              ></div>`
              : null
          }
        </div>
        ${
          hit.img_url
            ? html`<img
              src="${hit.img_url}"
              alt="${hit.name} spectrum"
              class="aa-ItemSpectra"
              loading="lazy"
            />`
            : null
        }
        ${hit.ex && hit.em ? html`<span class="aa-ItemInfo">${hit.ex}/${hit.em}</span>` : null}
      </div>
    </a>
  `
}

/**
 * Renders a reference search result
 * @param {Object} params
 * @param {Object} params.hit - The reference hit data
 * @param {Function} params.html - Algolia's HTML template function
 */
function ReferenceHit({ hit, html }) {
  const imageDir = IMAGE_DIR()

  return html`
    <a href="${hit.url}" class="aa-ItemLink">
      <article class="aa-ItemContent aa-ItemContent--reference">
        <div class="aa-ItemIcon">
          <img src="${imageDir}ref.png" alt="Reference icon" width="40" height="40" loading="lazy" />
        </div>
        <div class="aa-ItemContentBody">
          <div
            class="aa-ItemContentTitle"
            dangerouslySetInnerHTML=${{ __html: getHighlightedValue(hit, "citation") }}
          ></div>
          ${
            hit.title
              ? html`<div
                class="aa-ItemContentDescription"
                dangerouslySetInnerHTML=${{ __html: getHighlightedValue(hit, "title") }}
              ></div>`
              : null
          }
        </div>
      </article>
    </a>
  `
}

/**
 * Renders an organism search result
 * @param {Object} params
 * @param {Object} params.hit - The organism hit data
 * @param {Function} params.html - Algolia's HTML template function
 */
function OrganismHit({ hit, html }) {
  const imageDir = IMAGE_DIR()

  return html`
    <a href="${hit.url}" class="aa-ItemLink">
      <article class="aa-ItemContent aa-ItemContent--organism">
        <div class="aa-ItemIcon">
          <img
            src="${imageDir}organism_icon.png"
            alt="Organism icon"
            width="40"
            height="40"
            loading="lazy"
          />
        </div>
        <div class="aa-ItemContentBody">
          <em dangerouslySetInnerHTML=${{ __html: getHighlightedValue(hit, "scientific_name") }}></em>
        </div>
      </article>
    </a>
  `
}

// ============================================================================
// Source Factory Functions
// ============================================================================

/**
 * Creates a search source configuration
 * @param {Object} params
 * @param {string} params.sourceId - Unique identifier for the source
 * @param {string} params.indexName - Algolia index name
 * @param {string} params.query - Search query
 * @param {Object} params.searchClient - Algolia search client
 * @param {Object} params.config - Source-specific configuration
 * @param {Function} params.itemTemplate - Template function for rendering items
 * @param {string} params.headerLabel - Header label for the source section
 */
function createSearchSource({
  sourceId,
  indexName,
  query,
  searchClient,
  config,
  itemTemplate,
  headerLabel,
}) {
  return {
    sourceId,
    async getItems() {
      try {
        const { results } = await searchClient.search([
          {
            indexName,
            query,
            params: {
              hitsPerPage: config.hitsPerPage,
              attributesToRetrieve: config.attributes,
              attributesToHighlight: config.highlightAttributes,
            },
          },
        ])
        return results[0].hits
      } catch (error) {
        console.error(`Error searching ${sourceId}:`, error)
        return []
      }
    },
    templates: {
      header({ html }) {
        return html`<div class="aa-SourceHeader">${headerLabel}</div>`
      },
      item({ item, html }) {
        return itemTemplate({ hit: item, html })
      },
      ...(sourceId === "proteins" && {
        noResults({ html }) {
          return html`<div class="aa-ItemContent">No proteins found</div>`
        },
      }),
    },
    getItemUrl({ item }) {
      return item.url
    },
  }
}

/**
 * Creates all search sources
 * @param {string} query - The search query
 * @param {Object} searchClient - Algolia search client
 * @param {Object} algoliaConfig - Algolia configuration
 * @returns {Array} Array of search source configurations
 */
function createSearchSources(query, searchClient, algoliaConfig) {
  return [
    createSearchSource({
      sourceId: "proteins",
      indexName: algoliaConfig.proteinIndex,
      query,
      searchClient,
      config: SOURCE_CONFIG.proteins,
      itemTemplate: ProteinHit,
      headerLabel: "Proteins",
    }),
    createSearchSource({
      sourceId: "references",
      indexName: algoliaConfig.referenceIndex,
      query,
      searchClient,
      config: SOURCE_CONFIG.references,
      itemTemplate: ReferenceHit,
      headerLabel: "References",
    }),
    createSearchSource({
      sourceId: "organisms",
      indexName: algoliaConfig.organismIndex,
      query,
      searchClient,
      config: SOURCE_CONFIG.organisms,
      itemTemplate: OrganismHit,
      headerLabel: "Organisms",
    }),
  ]
}

// ============================================================================
// Plugin Factory Functions
// ============================================================================

/**
 * Creates the recent searches plugin
 * @returns {Object} Configured recent searches plugin
 */
function createRecentSearches() {
  return createLocalStorageRecentSearchesPlugin({
    key: SEARCH_CONFIG.recentSearchesKey,
    limit: SEARCH_CONFIG.recentSearchesLimit,
    transformSource({ source }) {
      return {
        ...source,
        templates: {
          ...source.templates,
          header({ html }) {
            return html`<div class="aa-SourceHeader">Recent Searches</div>`
          },
        },
      }
    },
  })
}

// ============================================================================
// Debug Helpers
// ============================================================================

/**
 * Activates debug mode by focusing input and triggering search
 * @param {HTMLElement} container - The autocomplete container
 */
function activateDebugMode(container) {
  setTimeout(() => {
    const input = container.querySelector(".aa-Input")
    if (input) {
      input.focus()
      input.dispatchEvent(new Event("input", { bubbles: true }))
    }
  }, 100)
}

// ============================================================================
// Main Component
// ============================================================================

/**
 * SearchAutocomplete - Algolia-powered search component
 *
 * @param {Object} props
 * @param {HTMLElement} props.container - The DOM container for the autocomplete
 * @param {boolean} [props.debug=false] - Enable debug mode (keeps panel open with test query)
 * @returns {null}
 */
export function SearchAutocomplete({ container, debug = false }) {
  const autocompleteRef = useRef(null)

  useEffect(() => {
    if (!container) return

    const algoliaConfig = window.FPBASE?.ALGOLIA
    if (!algoliaConfig) {
      console.error("Algolia configuration not found on window.FPBASE.ALGOLIA")
      return
    }

    const searchClient = liteClient(algoliaConfig.appID, algoliaConfig.publicKey)
    const recentSearchesPlugin = createRecentSearches()

    autocompleteRef.current = autocomplete({
      container,
      placeholder: SEARCH_CONFIG.placeholder,
      openOnFocus: true,
      plugins: [recentSearchesPlugin],
      detachedMediaQuery: SEARCH_CONFIG.detachedMediaQuery,
      insights: true,
      initialState: debug
        ? {
            query: SEARCH_CONFIG.debugQuery,
            isOpen: true,
          }
        : undefined,

      getSources({ query }) {
        if (!query && !debug) return []
        const searchQuery = query || SEARCH_CONFIG.debugQuery
        return createSearchSources(searchQuery, searchClient, algoliaConfig)
      },

      navigator: {
        navigate({ itemUrl }) {
          window.location.assign(itemUrl)
        },
        navigateNewTab({ itemUrl }) {
          const windowReference = window.open(itemUrl, "_blank", "noopener")
          windowReference?.focus()
        },
        navigateNewWindow({ itemUrl }) {
          window.open(itemUrl, "_blank", "noopener")
        },
      },
    })

    if (debug) {
      activateDebugMode(container)
    }

    return () => {
      autocompleteRef.current?.destroy()
    }
  }, [container, debug])

  return null
}

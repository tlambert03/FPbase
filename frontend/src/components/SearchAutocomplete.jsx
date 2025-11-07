import { autocomplete } from "@algolia/autocomplete-js"
import { createLocalStorageRecentSearchesPlugin } from "@algolia/autocomplete-plugin-recent-searches"
import { liteClient } from "algoliasearch/lite"
import { useEffect, useRef } from "react"
import "@algolia/autocomplete-theme-classic"

// Inject autocomplete styles
const styles = `
  /* Container */
  .autocomplete-container {
    width: 100%;
    flex: 1;
  }

  /* Input styling */
  .algolia-searchbar,
  .algolia-searchbar:focus {
    border: solid 1px rgba(255, 255, 255, 0.5);
    box-shadow: 0 1px 10px rgba(0, 0, 0, 0.2), 0 2px 4px 0 rgba(0, 0, 0, 0.1);
    color: #60a263;
    font-family: "Raleway", "Helvetica Neue", helvetica;
    height: 50px;
    padding-left: 16px;
    padding-right: 16px;
    border-radius: 4px;
    font-weight: 500;
    font-size: 1rem;
  }

  .algolia-searchbar::placeholder {
    color: #888;
    font-size: 1rem;
  }

  /* Autocomplete panel */
  .aa-Autocomplete,
  .aa-Form {
    width: 100%;
  }

  .aa-Panel {
    width: 600px !important;
    max-width: 90vw;
  }

  @media (max-width: 768px) {
    .aa-Panel {
      width: 100% !important;
    }
  }

  .aa-PanelLayout {
    max-height: 60vh;
    overflow-y: auto;
  }

  /* Source headers (Proteins, References, etc.) */
  .aa-SourceHeader {
    font-size: 0.85rem;
    font-weight: 600;
    color: #666;
    padding: 8px 12px 4px;
    border-top: 1px solid #eee;
  }

  .aa-SourceHeader:first-child {
    border-top: none;
  }

  /* Items */
  .aa-Item {
    border-bottom: 1px solid #f5f5f5;
  }

  .aa-Item:last-child {
    border-bottom: none;
  }

  .aa-Item[aria-selected="true"] .aa-ItemContent {
    background-color: #eee;
  }

  .aa-ItemLink {
    text-decoration: none;
    color: inherit;
    display: block;
  }

  /* Item content grid layout */
  .aa-ItemContent {
    display: grid;
    grid-template-columns: 32px 1fr auto auto;
    align-items: center;
    gap: 8px;
    padding: 6px 12px;
    min-height: 44px;
  }

  .aa-ItemContent:hover {
    background-color: #f5f5f5;
  }

  .aa-ItemIcon {
    grid-column: 1;
    width: 32px;
    height: 32px;
    opacity: 0.75;
  }

  .aa-ItemIcon--protein {
    opacity: 0.9;
  }

  .aa-ItemBody {
    grid-column: 2;
    min-width: 0;
    line-height: 1.4;
  }

  .aa-ItemTitle {
    font-weight: 500;
    font-size: 0.95rem;
  }

  .aa-ItemTitle em {
    font-weight: 700;
    font-style: normal;
  }

  .aa-ItemAliases {
    font-size: 0.75rem;
    color: #777;
    margin-top: 1px;
  }

  .aa-ItemAliases em {
    font-weight: 600;
    font-style: normal;
  }

  .aa-ItemInfo {
    grid-column: 3;
    font-size: 0.7rem;
    color: #888;
    white-space: nowrap;
    padding-right: 8px;
  }

  .aa-ItemSpectra {
    grid-column: 4;
    width: 120px;
    height: 36px;
    object-fit: contain;
    opacity: 0.7;
    filter: grayscale(50%);
  }

  .aa-Item[aria-selected="true"] .aa-ItemSpectra {
    opacity: 0.95;
    filter: grayscale(20%);
  }
`

if (typeof document !== "undefined") {
  const styleId = "autocomplete-styles"
  if (!document.getElementById(styleId)) {
    const styleTag = document.createElement("style")
    styleTag.id = styleId
    styleTag.textContent = styles
    document.head.appendChild(styleTag)
  }
}

function ProteinHit({ hit, html }) {
  const imageDir = window.FPBASE?.imageDir || "/static/images/"

  let iconColor = "gray50"
  if (hit.switchType && hit.switchType !== "Basic") {
    iconColor = "rainbow"
  } else if (hit.color && !hit.color.includes("Stokes")) {
    iconColor = hit.color.toLowerCase().replace(/ |\//g, "_")
  }

  return html`
    <a href="${hit.url}" class="aa-ItemLink">
      <div class="aa-ItemContent">
        <img
          src="${imageDir}gfp_${iconColor}_40.png"
          alt="Protein icon"
          class="aa-ItemIcon aa-ItemIcon--protein"
          width="40"
          height="40"
          loading="lazy"
        />
        <div class="aa-ItemBody">
          <div
            class="aa-ItemTitle"
            dangerouslySetInnerHTML=${{
              __html: hit._highlightResult?.name?.value || hit.name,
            }}
          ></div>
          ${
            hit.aliases && hit.aliases.length > 0
              ? html`<div
                  class="aa-ItemAliases"
                  dangerouslySetInnerHTML=${{
                    __html:
                      "aka: " +
                      (hit._highlightResult?.aliases
                        ? hit._highlightResult.aliases
                            .slice(0, 3)
                            .map((a) => a.value)
                            .join(", ")
                        : hit.aliases.slice(0, 3).join(", ")),
                  }}
                ></div>`
              : null
          }
        </div>
        ${hit.ex && hit.em ? html`<span class="aa-ItemInfo">${hit.ex}/${hit.em}</span>` : null}
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
      </div>
    </a>
  `
}

function ReferenceHit({ hit, html }) {
  const imageDir = window.FPBASE?.imageDir || "/static/images/"

  return html`
    <a href="${hit.url}" class="aa-ItemLink">
      <article class="aa-ItemContent">
        <div class="aa-ItemIcon">
          <img
            src="${imageDir}ref.png"
            alt="Reference icon"
            width="40"
            height="40"
            loading="lazy"
          />
        </div>
        <div class="aa-ItemContentBody">
          <div
            class="aa-ItemContentTitle"
            dangerouslySetInnerHTML=${{
              __html: hit._highlightResult?.citation?.value || hit.citation,
            }}
          ></div>
          ${
            hit.title
              ? html`
            <div
              class="aa-ItemContentDescription"
              dangerouslySetInnerHTML=${{
                __html: hit._highlightResult?.title?.value || hit.title,
              }}
            ></div>
          `
              : null
          }
        </div>
      </article>
    </a>
  `
}

function OrganismHit({ hit, html }) {
  const imageDir = window.FPBASE?.imageDir || "/static/images/"

  return html`
    <a href="${hit.url}" class="aa-ItemLink">
      <article class="aa-ItemContent">
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
          <em
            dangerouslySetInnerHTML=${{
              __html: hit._highlightResult?.scientific_name?.value || hit.scientific_name,
            }}
          ></em>
        </div>
      </article>
    </a>
  `
}

export function SearchAutocomplete({ container }) {
  const autocompleteRef = useRef(null)

  useEffect(() => {
    if (!container) return

    // Initialize Algolia config and clients inside useEffect to avoid module-level window access
    const ALGOLIA_CONFIG = window.FPBASE?.ALGOLIA
    if (!ALGOLIA_CONFIG) {
      console.error("Algolia configuration not found on window.FPBASE.ALGOLIA")
      return
    }

    const searchClient = liteClient(ALGOLIA_CONFIG.appID, ALGOLIA_CONFIG.publicKey)

    // Recent searches plugin
    const recentSearchesPlugin = createLocalStorageRecentSearchesPlugin({
      key: "fpbase-search",
      limit: 5,
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

    autocompleteRef.current = autocomplete({
      container,
      placeholder: "Search proteins, organisms, references...",
      openOnFocus: true,
      plugins: [recentSearchesPlugin],
      detachedMediaQuery: "(max-width: 768px)", // Mobile: full-screen overlay
      insights: true, // Enable insights for analytics

      getSources({ query }) {
        if (!query) return []

        return [
          // Proteins
          {
            sourceId: "proteins",
            getItems() {
              return searchClient
                .search([
                  {
                    indexName: ALGOLIA_CONFIG.proteinIndex,
                    query,
                    params: {
                      hitsPerPage: 5,
                      attributesToRetrieve: [
                        "name",
                        "uuid",
                        "aliases",
                        "switchType",
                        "color",
                        "ex",
                        "em",
                        "img_url",
                        "url",
                      ],
                      attributesToHighlight: ["name", "aliases"],
                    },
                  },
                ])
                .then(({ results }) => results[0].hits)
                .catch((error) => {
                  console.error("Error searching proteins:", error)
                  return []
                })
            },
            templates: {
              header({ html }) {
                return html`<div class="aa-SourceHeader">Proteins</div>`
              },
              item({ item, html }) {
                return ProteinHit({ hit: item, html })
              },
              noResults({ html }) {
                return html`<div class="aa-ItemContent">No proteins found</div>`
              },
            },
            getItemUrl({ item }) {
              return item.url
            },
          },

          // References
          {
            sourceId: "references",
            getItems() {
              return searchClient
                .search([
                  {
                    indexName: ALGOLIA_CONFIG.referenceIndex,
                    query,
                    params: {
                      hitsPerPage: 3,
                      attributesToHighlight: ["citation", "title"],
                    },
                  },
                ])
                .then(({ results }) => results[0].hits)
                .catch((error) => {
                  console.error("Error searching references:", error)
                  return []
                })
            },
            templates: {
              header({ html }) {
                return html`<div class="aa-SourceHeader">References</div>`
              },
              item({ item, html }) {
                return ReferenceHit({ hit: item, html })
              },
            },
            getItemUrl({ item }) {
              return item.url
            },
          },

          // Organisms
          {
            sourceId: "organisms",
            getItems() {
              return searchClient
                .search([
                  {
                    indexName: ALGOLIA_CONFIG.organismIndex,
                    query,
                    params: {
                      hitsPerPage: 2,
                      attributesToHighlight: ["scientific_name"],
                    },
                  },
                ])
                .then(({ results }) => results[0].hits)
                .catch((error) => {
                  console.error("Error searching organisms:", error)
                  return []
                })
            },
            templates: {
              header({ html }) {
                return html`<div class="aa-SourceHeader">Organisms</div>`
              },
              item({ item, html }) {
                return OrganismHit({ hit: item, html })
              },
            },
            getItemUrl({ item }) {
              return item.url
            },
          },
        ]
      },

      navigator: {
        navigate({ itemUrl }) {
          window.location.assign(itemUrl)
        },
        navigateNewTab({ itemUrl }) {
          const windowReference = window.open(itemUrl, "_blank", "noopener")
          if (windowReference) {
            windowReference.focus()
          }
        },
        navigateNewWindow({ itemUrl }) {
          window.open(itemUrl, "_blank", "noopener")
        },
      },
    })

    return () => {
      autocompleteRef.current?.destroy()
    }
  }, [container])

  return null
}

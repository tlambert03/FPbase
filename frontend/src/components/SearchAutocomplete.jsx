import { autocomplete } from "@algolia/autocomplete-js"
import { createLocalStorageRecentSearchesPlugin } from "@algolia/autocomplete-plugin-recent-searches"
import { liteClient } from "algoliasearch/lite"
import { useEffect, useRef } from "react"
import "@algolia/autocomplete-theme-classic"

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
      <article class="aa-ItemContent">
        <div class="aa-ItemIcon aa-ItemIcon--protein">
          <img
            src="${imageDir}gfp_${iconColor}_40.png"
            alt="Protein icon"
            width="40"
            height="40"
            loading="lazy"
          />
        </div>
        <div class="aa-ItemContentBody">
          <div class="aa-ItemContentTitle">
            ${hit._highlightResult?.name?.value || hit.name}
          </div>
          ${
            hit.aliases && hit.aliases.length > 0
              ? html`
            <div class="aa-ItemContentDescription">
              aka: ${
                hit._highlightResult?.aliases
                  ? hit._highlightResult.aliases
                      .slice(0, 3)
                      .map((a) => a.value)
                      .join(", ")
                  : hit.aliases.slice(0, 3).join(", ")
              }
            </div>
          `
              : null
          }
          ${
            hit.ex && hit.em
              ? html`
            <div class="aa-ItemContentSubtitle">
              <abbr title="Excitation/Emission wavelengths">${hit.ex}/${hit.em} nm</abbr>
            </div>
          `
              : null
          }
        </div>
        ${
          hit.img_url
            ? html`
          <img
            src="${hit.img_url}"
            alt="${hit.name} spectrum"
            class="aa-ItemContentSpectra"
            loading="lazy"
          />
        `
            : null
        }
      </article>
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
          <div class="aa-ItemContentTitle">
            ${hit._highlightResult?.citation?.value || hit.citation}
          </div>
          ${
            hit.title
              ? html`
            <div class="aa-ItemContentDescription">
              ${hit._highlightResult?.title?.value || hit.title}
            </div>
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
          <em>${hit._highlightResult?.scientific_name?.value || hit.scientific_name}</em>
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

import { useEffect, useRef } from 'react';
import { autocomplete } from '@algolia/autocomplete-js';
import { createLocalStorageRecentSearchesPlugin } from '@algolia/autocomplete-plugin-recent-searches';
import { liteClient as algoliasearch } from 'algoliasearch/lite';
import '@algolia/autocomplete-theme-classic';

function ProteinHit({ hit, html }) {
  const imageDir = window.FPBASE?.imageDir || '/static/images/';

  let iconColor = 'gray50';
  if (hit.switchType && hit.switchType !== 'Basic') {
    iconColor = 'rainbow';
  } else if (hit.color && !hit.color.includes('Stokes')) {
    iconColor = hit.color.toLowerCase().replace(/ |\//g, '_');
  }

  return html`
    <a href="${hit.url}" class="aa-ItemLink">
      <div class="aa-ItemContent">
        <div class="aa-ItemIcon aa-ItemIcon--protein">
          <img src="${imageDir}gfp_${iconColor}_40.png" alt="" />
        </div>
        <div class="aa-ItemContentBody">
          <div class="aa-ItemContentTitle">${hit.name}</div>
          ${hit.aliases && hit.aliases.length > 0 ? html`
            <div class="aa-ItemContentDescription">aka: ${hit.aliases.slice(0, 3).join(', ')}</div>
          ` : null}
          ${hit.ex && hit.em ? html`
            <div class="aa-ItemContentSubtitle">${hit.ex}/${hit.em} nm</div>
          ` : null}
        </div>
        ${hit.img_url ? html`
          <img src="${hit.img_url}" alt="spectrum" class="aa-ItemContentSpectra" />
        ` : null}
      </div>
    </a>
  `;
}

function ReferenceHit({ hit, html }) {
  const imageDir = window.FPBASE?.imageDir || '/static/images/';

  return html`
    <a href="${hit.url}" class="aa-ItemLink">
      <div class="aa-ItemContent">
        <div class="aa-ItemIcon">
          <img src="${imageDir}ref.png" alt="" />
        </div>
        <div class="aa-ItemContentBody">
          <div class="aa-ItemContentTitle">${hit.citation}</div>
          ${hit.title ? html`
            <div class="aa-ItemContentDescription">${hit.title}</div>
          ` : null}
        </div>
      </div>
    </a>
  `;
}

function OrganismHit({ hit, html }) {
  const imageDir = window.FPBASE?.imageDir || '/static/images/';

  return html`
    <a href="${hit.url}" class="aa-ItemLink">
      <div class="aa-ItemContent">
        <div class="aa-ItemIcon">
          <img src="${imageDir}organism_icon.png" alt="" />
        </div>
        <div class="aa-ItemContentBody">
          ${hit.scientific_name}
        </div>
      </div>
    </a>
  `;
}

export function SearchAutocomplete({ container }) {
  const autocompleteRef = useRef(null);

  useEffect(() => {
    if (!container) return;

    // Initialize Algolia config and clients inside useEffect to avoid module-level window access
    const ALGOLIA_CONFIG = window.FPBASE?.ALGOLIA;
    if (!ALGOLIA_CONFIG) {
      console.error('Algolia configuration not found on window.FPBASE.ALGOLIA');
      return;
    }

    const searchClient = algoliasearch(ALGOLIA_CONFIG.appID, ALGOLIA_CONFIG.publicKey);

    // Recent searches plugin
    const recentSearchesPlugin = createLocalStorageRecentSearchesPlugin({
      key: 'fpbase-search',
      limit: 5,
      transformSource({ source }) {
        return {
          ...source,
          templates: {
            ...source.templates,
            header({ html }) {
              return html`<div class="aa-SourceHeader">Recent Searches</div>`;
            },
          },
        };
      },
    });

    autocompleteRef.current = autocomplete({
      container,
      placeholder: 'Search proteins, organisms, references...',
      openOnFocus: true,
      plugins: [recentSearchesPlugin],
      detachedMediaQuery: '(max-width: 768px)', // Mobile: full-screen overlay

      getSources({ query }) {
        if (!query) return [];

        return [
          // Proteins
          {
            sourceId: 'proteins',
            getItems() {
              return searchClient.search([{
                indexName: ALGOLIA_CONFIG.proteinIndex,
                query,
                params: {
                  hitsPerPage: 5,
                  attributesToRetrieve: [
                    'name', 'uuid', 'aliases', 'switchType', 'color',
                    'ex', 'em', 'img_url', 'url'
                  ],
                },
              }]).then(({ results }) => results[0].hits);
            },
            templates: {
              header({ html }) {
                return html`<div class="aa-SourceHeader">Proteins</div>`;
              },
              item({ item, html }) {
                return ProteinHit({ hit: item, html });
              },
              noResults({ html }) {
                return html`<div class="aa-ItemContent">No proteins found</div>`;
              },
            },
            getItemUrl({ item }) {
              return item.url;
            },
          },

          // References
          {
            sourceId: 'references',
            getItems() {
              return searchClient.search([{
                indexName: ALGOLIA_CONFIG.referenceIndex,
                query,
                params: { hitsPerPage: 3 },
              }]).then(({ results }) => results[0].hits);
            },
            templates: {
              header({ html }) {
                return html`<div class="aa-SourceHeader">References</div>`;
              },
              item({ item, html }) {
                return ReferenceHit({ hit: item, html });
              },
            },
            getItemUrl({ item }) {
              return item.url;
            },
          },

          // Organisms
          {
            sourceId: 'organisms',
            getItems() {
              return searchClient.search([{
                indexName: ALGOLIA_CONFIG.organismIndex,
                query,
                params: { hitsPerPage: 2 },
              }]).then(({ results }) => results[0].hits);
            },
            templates: {
              header({ html }) {
                return html`<div class="aa-SourceHeader">Organisms</div>`;
              },
              item({ item, html }) {
                return OrganismHit({ hit: item, html });
              },
            },
            getItemUrl({ item }) {
              return item.url;
            },
          },
        ];
      },

      navigator: {
        navigate({ itemUrl }) {
          window.location.assign(itemUrl);
        },
      },
    });

    return () => {
      autocompleteRef.current?.destroy();
    };
  }, [container]);

  return null;
}

import { useEffect, useRef } from 'react';
import { autocomplete } from '@algolia/autocomplete-js';
import { createLocalStorageRecentSearchesPlugin } from '@algolia/autocomplete-plugin-recent-searches';
import { liteClient as algoliasearch } from 'algoliasearch/lite';
import '@algolia/autocomplete-theme-classic';

function ProteinHit({ hit, components }) {
  const imageDir = window.FPBASE?.imageDir || '/static/images/';

  let iconColor = 'gray50';
  if (hit.switchType && hit.switchType !== 'Basic') {
    iconColor = 'rainbow';
  } else if (hit.color && !hit.color.includes('Stokes')) {
    iconColor = hit.color.toLowerCase().replace(/ |\//g, '_');
  }

  return (
    <a href={hit.url} className="aa-ItemLink">
      <div className="aa-ItemContent">
        <div className="aa-ItemIcon aa-ItemIcon--protein">
          <img
            src={`${imageDir}gfp_${iconColor}_40.png`}
            alt=""
          />
        </div>
        <div className="aa-ItemContentBody">
          <div className="aa-ItemContentTitle">
            <components.Highlight hit={hit} attribute="name" />
          </div>
          {hit.aliases && hit.aliases.length > 0 && (
            <div className="aa-ItemContentDescription">
              aka: {hit.aliases.slice(0, 3).join(', ')}
            </div>
          )}
          {hit.ex && hit.em && (
            <div className="aa-ItemContentSubtitle">
              {hit.ex}/{hit.em} nm
            </div>
          )}
        </div>
        {hit.img_url && (
          <img
            src={hit.img_url}
            alt="spectrum"
            className="aa-ItemContentSpectra"
          />
        )}
      </div>
    </a>
  );
}

function ReferenceHit({ hit, components }) {
  const imageDir = window.FPBASE?.imageDir || '/static/images/';

  return (
    <a href={hit.url} className="aa-ItemLink">
      <div className="aa-ItemContent">
        <div className="aa-ItemIcon">
          <img src={`${imageDir}ref.png`} alt="" />
        </div>
        <div className="aa-ItemContentBody">
          <div className="aa-ItemContentTitle">
            <components.Highlight hit={hit} attribute="citation" />
          </div>
          {hit._highlightResult?.title && (
            <div className="aa-ItemContentDescription">
              <components.Highlight hit={hit} attribute="title" />
            </div>
          )}
        </div>
      </div>
    </a>
  );
}

function OrganismHit({ hit, components }) {
  const imageDir = window.FPBASE?.imageDir || '/static/images/';

  return (
    <a href={hit.url} className="aa-ItemLink">
      <div className="aa-ItemContent">
        <div className="aa-ItemIcon">
          <img src={`${imageDir}organism_icon.png`} alt="" />
        </div>
        <div className="aa-ItemContentBody">
          <components.Highlight hit={hit} attribute="scientific_name" />
        </div>
      </div>
    </a>
  );
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
            header() {
              return <div className="aa-SourceHeader">Recent Searches</div>;
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
              header() {
                return <div className="aa-SourceHeader">Proteins</div>;
              },
              item({ item, components }) {
                return <ProteinHit hit={item} components={components} />;
              },
              noResults() {
                return <div className="aa-ItemContent">No proteins found</div>;
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
              header() {
                return <div className="aa-SourceHeader">References</div>;
              },
              item({ item, components }) {
                return <ReferenceHit hit={item} components={components} />;
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
              header() {
                return <div className="aa-SourceHeader">Organisms</div>;
              },
              item({ item, components }) {
                return <OrganismHit hit={item} components={components} />;
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

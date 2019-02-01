function checkObject(val, prop, str){
  var propDict = {
    'genbank': 'GenBank',
    'pdb': 'PDB',
    'ipg_id': 'IPG ID',
    'uniprot': 'UniProt',
    'aliases': 'aka',
    'doi': 'DOI',
    'pmid': 'PMID',
    'title': 'Title',
    '_excerpts': 'Excerpt',
    'prot_primary': 'Protein',
    'prot_secondary': 'Protein'
  }
  if (val.matchLevel === 'full' || val.matchLevel === 'partial'){
    if (str.length > 0) { str = str + "; "}
    str = str + propDict[prop] + ": " + val.value;
  }
  return str
}

function highlightHits(high) {
  var str = '';
  for (var prop in high) {
    if (high.hasOwnProperty(prop) && prop !== 'name' && prop !== 'citation' ) {
      if (high[prop].constructor === Array) {
        for (var i = 0; i < high[prop].length; i++) {
          if (typeof high[prop][i] === 'object'){
            str = checkObject(high[prop][i], prop, str)
          }
        }
      } else {
        if (typeof high[prop] === 'object'){
          str = checkObject(high[prop], prop, str)
        }
      }
    }
  }
  if (str){
    return "<span class='highlighted-hits'>(" + str + ")</span>";
  } else {
    return ''
  }
}

function highlightRefHits(high) {
  var str = '';
  for (var prop in high) {
    if (high.hasOwnProperty(prop) && prop !== 'name' && prop !== 'citation' ) {
      if (high[prop].constructor === Array) {
        for (var i = 0; i < high[prop].length; i++) {
          if (typeof high[prop][i] === 'object'){
            str = checkObject(high[prop][i], prop, str)
          }
        }
      } else {
        if (typeof high[prop] === 'object'){
          str = checkObject(high[prop], prop, str)
        }
      }
    }
  }
  if (str){
    return "<span class='highlighted-hits'>(" + str + ")</span>";
  } else {
    return ''
  }
}


export default function initAutocomplete() {
  var algoliaClient = algoliasearch(FPBASE.ALGOLIA.appID, FPBASE.ALGOLIA.publicKey);
  var proteinIndex = algoliaClient.initIndex(FPBASE.ALGOLIA.proteinIndex);
  var organismIndex = algoliaClient.initIndex(FPBASE.ALGOLIA.organismIndex);
  var referenceIndex = algoliaClient.initIndex(FPBASE.ALGOLIA.referenceIndex);
  autocomplete('#algolia-search-input',
    {
      getRankingInfo: false,
      minLength: 3,
      autoselect: false,
      templates: {
        footer: '<div class="branding">search by <a href="https://algolia.com"><img src="https://www.algolia.com/static_assets/images/press/downloads/algolia-logo-light.svg" /></a></div>'
      }
    },
    [
      {
        // add {attributesToRetrieve: }
        source: autocomplete.sources.hits(proteinIndex, { hitsPerPage: 6, getRankingInfo: 1}),
        displayKey: 'name',
        templates: {
          suggestion: function(suggestion) {
            if (suggestion.switchType && suggestion.switchType !== 'Basic') {
              var col = 'rainbow';
            } else if (suggestion.color && !suggestion.color.includes('Stokes') ) {
              var col = suggestion.color.toLowerCase().replace(/ |\//g,"_");
            } else {
              var col = 'gray50'
            }
            var str = "<img class='type protein' src='" + FPBASE.imageDir + "gfp_" + col + "_40.png'>";
            str = str + suggestion._highlightResult.name.value;
            if (suggestion.img_url){
              str = str + "<img class='spectra' src=" + suggestion.img_url + ">";
            }
            str = str + highlightHits(suggestion._highlightResult)
            var info = ''
            if (suggestion.switchType === 'Basic') {
              if (suggestion.ex && suggestion.em){
                info = suggestion.ex + "/" + suggestion.em
              } else {
                info = ''
              }
            } else if (suggestion.switchType) {
              info = {
                'photoswitchable': 'PS',
                'photoactivatable': 'PA',
                'photoconvertible': 'PC',
                'multi-photochromic': 'MPC',
                'multistate': 'MS',
                'timer': 'Time'
              }[suggestion.switchType.toLowerCase()]
            }
            str = str + "<span class='info'>" + info + "</span>";
            return "<a href='" + suggestion.url + "'><div>" + str + "</div></a>";
          }
        }
      },
      {
        source: autocomplete.sources.hits(referenceIndex, { hitsPerPage: 2, getRankingInfo: 1}),
        displayKey: 'citation',
        templates: {
          suggestion: function(suggestion) {
            var str = suggestion._highlightResult.citation.value;
            str = str + "<img class='type' src='" + FPBASE.imageDir + "ref.png" + "'>";
            str = str + highlightRefHits(suggestion._highlightResult)
            return "<a href='" + suggestion.url + "'><div>" + str + "</div></a>";
          }
        }
      },
      {
        source: autocomplete.sources.hits(organismIndex, { hitsPerPage: 2, getRankingInfo: 1 }),
        displayKey: 'scientific_name',
        templates: {
          suggestion: function(suggestion) {
            var str = suggestion._highlightResult.scientific_name.value;
            str = str + "<img class='type' src='" + FPBASE.imageDir + "organism_icon.png" + "'>";
            return "<a href='" + suggestion.url + "'><div>" + str + "</div></a>";
          }
        }
      }
    ]
  ).on('autocomplete:selected', function(event, suggestion, dataset, context) {
    console.log(suggestion, dataset);
    if (context.selectionMethod === 'click') {
      return;
    }
    // Change the page, for example, on other events
    window.location.assign(suggestion.url);
  }).on('autocomplete:change', function(event, suggestion, dataset, context) {
    console.log('changed')
  })
}

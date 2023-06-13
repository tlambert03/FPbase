import $ from "jquery"

function checkObject(val, prop, str) {
  var propDict = {
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
    prot_secondary: "Protein"
  }
  if (val.matchLevel === "full" || val.matchLevel === "partial") {
    if (str.length > 0) {
      str = str + "; "
    }
    str = str + propDict[prop] + ": " + val.value
  }
  return str
}

function highlightHits(high) {
  var str = ""
  for (var prop in high) {
    if (high.hasOwnProperty(prop) && prop !== "name" && prop !== "citation") {
      if (high[prop].constructor === Array) {
        for (var i = 0; i < high[prop].length; i++) {
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
    return "<span class='highlighted-hits'>(" + str + ")</span>"
  } else {
    return ""
  }
}

// function highlightRefHits(high) {
//   var str = '';
//   for (var prop in high) {
//     if (high.hasOwnProperty(prop) && prop !== 'name' && prop !== 'citation' ) {
//       if (high[prop].constructor === Array) {
//         for (var i = 0; i < high[prop].length; i++) {
//           if (typeof high[prop][i] === 'object'){
//             str = checkObject(high[prop][i], prop, str)
//           }
//         }
//       } else {
//         if (typeof high[prop] === 'object'){
//           str = checkObject(high[prop], prop, str)
//         }
//       }
//     }
//   }
//   if (str){
//     return "<span class='highlighted-hits'>(" + str + ")</span>";
//   } else {
//     return ''
//   }
// }

function highlightRefHits(high) {
  function recurseMatches(obj) {
    var results = {}

    function innerRecurse(obj, _key) {
      if (_key !== undefined) {
        obj = obj[_key]
      }
      for (var key in obj) {
        if (obj.hasOwnProperty(key)) {
          if (obj[key].constructor === Array) {
            innerRecurse(obj, key)
          } else if (obj[key].hasOwnProperty("matchLevel")) {
            if (
              obj[key].matchLevel !== "none" &&
              obj[key].value !== undefined
            ) {
              let attr = key
              if (_key !== undefined) {
                attr = _key
              }
              if (!results.hasOwnProperty(attr)) {
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
    return (
      str
        .split(" ")
        .splice(0, no_words)
        .join(" ") + " ..."
    )
  }

  var results = recurseMatches(high)

  if (results.hasOwnProperty("prot_primary")) {
    delete results.prot_secondary
  }
  var str = ""
  var items = [
    ["doi", "DOI"],
    ["pmid", "PMID"],
    ["prot_primary", "Protein"],
    ["prot_secondary", "2˚ Protein"]
  ]
  for (var x = 0; x < items.length; x++) {
    var key = items[x][0]
    var title = items[x][1]
    if (results.hasOwnProperty(key)) {
      var _str = []
      for (var i = 0; i < results[key].length; i++) {
        _str.push(results[key][i][0])
      }
      str = str + title + ": " + _str.join(", ")
    }
  }
  if (str) {
    str = "<span class='highlighted-hits'>(" + str + ")</span>"
  }
  if (results.hasOwnProperty("title")) {
    str = str + '<div class="ref-title" >' + results.title[0][0] + "</div>"
  } else if (results.hasOwnProperty("_excerpts")) {
    if (
      results._excerpts.some(function(d) {
        return d[1] === "full"
      })
    ) {
      for (var e = 0; e < results._excerpts.length; e++) {
        if (results._excerpts[e][1] === "full") {
          var exc = results._excerpts[e][0]
          var pre = exc.split("<em>")[0].split(" ")
          var _pre = ""
          if (pre.length > 5) {
            _pre = _pre + "... "
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
    } else {
    }
  }
  return str
}

export default async function initAutocomplete() {
  const [{ default: algoliasearch }] = await Promise.all([
    import("algoliasearch"),
    import("autocomplete.js/dist/autocomplete.jquery.js")
  ])

  var algoliaClient = algoliasearch(
    window.FPBASE.ALGOLIA.appID,
    window.FPBASE.ALGOLIA.publicKey
  )
  var proteinIndex = algoliaClient.initIndex(window.FPBASE.ALGOLIA.proteinIndex)
  var organismIndex = algoliaClient.initIndex(
    window.FPBASE.ALGOLIA.organismIndex
  )
  var referenceIndex = algoliaClient.initIndex(
    window.FPBASE.ALGOLIA.referenceIndex
  )

  function empty(context) {
    if (context.hasOwnProperty("query")) {
      var p = context.query.trim().replace(/\s/g, "%20")
      return (
        '<div class="empty"><span class="nohits"></span>No results... try the <a href="/search/?name__icontains=' +
        p +
        '">advanced search</a></div>'
      )
    } else {
      return '<div class="empty"><span class="nohits"></span>No results... try the <a href="/search/">advanced search</a></div>'
    }
  }
  $("#algolia-search-input")
    .autocomplete(
      {
        getRankingInfo: false,
        minLength: 3,
        autoselect: true,
        autoselectOnBlur: window.mobilecheck(),
        templates: {
          empty: empty
        }
      },
      [
        {
          // add {attributesToRetrieve: }
          source: $.fn.autocomplete.sources.hits(proteinIndex, {
            hitsPerPage: 5
          }),
          displayKey: "name",
          templates: {
            suggestion: function(suggestion) {
              let col
              if (suggestion.switchType && suggestion.switchType !== "Basic") {
                col = "rainbow"
              } else if (
                suggestion.color &&
                !suggestion.color.includes("Stokes")
              ) {
                col = suggestion.color.toLowerCase().replace(/ |\//g, "_")
              } else {
                col = "gray50"
              }
              var str =
                "<img class='type protein' src='" +
                window.FPBASE.imageDir +
                "gfp_" +
                col +
                "_40.png'>"
              str = str + suggestion._highlightResult.name.value
              if (suggestion.img_url) {
                str =
                  str + "<img class='spectra' src=" + suggestion.img_url + ">"
              }
              str = str + highlightHits(suggestion._highlightResult)
              var info = ""
              if (suggestion.switchType === "Basic") {
                if (suggestion.ex && suggestion.em) {
                  info = suggestion.ex + "/" + suggestion.em
                } else {
                  info = ""
                }
              } else if (suggestion.switchType) {
                info = {
                  photoswitchable: "PS",
                  photoactivatable: "PA",
                  photoconvertible: "PC",
                  "multi-photochromic": "MPC",
                  multistate: "MS",
                  timer: "Time"
                }[suggestion.switchType.toLowerCase()]
              }
              str = str + "<span class='info'>" + info + "</span>"
              return (
                "<a href='" + suggestion.url + "'><div>" + str + "</div></a>"
              )
            }
          }
        },
        {
          source: $.fn.autocomplete.sources.hits(referenceIndex, {
            hitsPerPage: 3
          }),
          displayKey: "citation",
          templates: {
            suggestion: function(suggestion) {
              var str = suggestion._highlightResult.citation.value
              str =
                str +
                "<img class='type' src='" +
                window.FPBASE.imageDir +
                "ref.png" +
                "'>"
              str = str + highlightRefHits(suggestion._highlightResult)
              return (
                "<a href='" + suggestion.url + "'><div>" + str + "</div></a>"
              )
            }
          }
        },
        {
          source: $.fn.autocomplete.sources.hits(organismIndex, {
            hitsPerPage: 2
          }),
          displayKey: "scientific_name",
          templates: {
            suggestion: function(suggestion) {
              var str = suggestion._highlightResult.scientific_name.value
              str =
                str +
                "<img class='type' src='" +
                window.FPBASE.imageDir +
                "organism_icon.png" +
                "'>"
              return (
                "<a href='" + suggestion.url + "'><div>" + str + "</div></a>"
              )
            }
          }
        },
        {
          source: function(query, callback) {
            callback([
              {
                query: query,
                url: "/search/?q=" + encodeURI(query)
              }
            ])
          },
          templates: {
            suggestion: function(suggestion) {
              return (
                '<div class="search-footer"><a class="asearch" href="' +
                suggestion.url +
                '">Advanced search for: <em>' +
                suggestion.query +
                '</em></a><div class="branding">search powered by <a href="https://algolia.com"><img src="' +
                window.FPBASE.imageDir +
                'logo-algolia-nebula-blue-full.svg" /></a></div></div>'
              )
            }
          }
        }
      ]
    )
    .on("autocomplete:selected", function(event, suggestion, dataset, context) {
      if (context.selectionMethod === "click") {
        return
      }
      // Change the page, for example, on other events
      window.location.assign(suggestion.url)
    })
}

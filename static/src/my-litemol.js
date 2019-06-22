import "./css/litemol/LiteMol-plugin-blue.css"
import LiteMol from "./js/pdb/LiteMol-plugin.js"
import $ from "jquery"

function initLiteMol(selection, changer) {
  const PluginSpec = LiteMol.Plugin.getDefaultSpecification()
  const LayoutRegion = LiteMol.Bootstrap.Components.LayoutRegion
  const Components = LiteMol.Plugin.Components
  PluginSpec.components = [
    Components.Visualization.HighlightInfo(LayoutRegion.Main, true),
    Components.Entity.Current("LiteMol", LiteMol.Plugin.VERSION.number)(
      LayoutRegion.Right,
      true
    ),
    Components.Transform.View(LayoutRegion.Right),
    //Components.Context.Log(LayoutRegion.Bottom, true),
    Components.Context.Overlay(LayoutRegion.Root),
    //Components.Context.Toast(LayoutRegion.Main, true),
    Components.Context.BackgroundTasks(LayoutRegion.Main, true)
  ]

  try {
    var plugin = LiteMol.Plugin.create({
      customSpecification: PluginSpec,
      target: selection,
      viewportBackground: "#fff",
      layoutState: {
        hideControls: true,
        isExpanded: false
      },
      allowAnalytics: true
    })

    var dataCache = {}
    changer.change(function() {
      var id = this.value
      if (!dataCache.hasOwnProperty(id)) {
        dataCache[id] = getPDBbinary(id)
      }
      plugin.clear()
      dataCache[id].then(
        data =>
          plugin.loadMolecule({
            data,
            id
          }),
        reason => {
          $(selection).html(
            '<span class="text-danger muted">failed to retrieve molecular structure</span>'
          )
        }
      )
    })

    $("body").click(function(e) {
      if ($(".lm-layout-right").length) {
        if ($(e.target).closest("#litemol-viewer").length === 0) {
          plugin.setLayoutState({
            hideControls: true
          })
        }
      }
    })
  } catch (err) {
    if (window.Sentry !== undefined) {
      window.Sentry.captureException(err)
    }
  }

  changer
    .change(function() {
      // var id = this.value
      $("#pdb-out-link").attr(
        "href",
        "https://www.rcsb.org/structure/" + this.value
      )
      if (pdb_info[this.value]) {
        loadSmiles(this.value)
        loadChemInfo(this.value)
      }
    })
    .trigger("change")
}

function getPDBbinary(id) {
  return new Promise(function(resolve, reject) {
    $.get(
      "https://www.ebi.ac.uk/pdbe/static/entry/" +
        id.toLowerCase() +
        "_updated.cif"
    )
      .done(response => {
        resolve(response)
      })
      .fail((xhr, status, error) => {
        $.get("https://files.rcsb.org/download/" + id + ".cif")
          .done(response => {
            resolve(response)
          })
          .fail(xhr =>
            reject({
              status: xhr.status,
              statusText: xhr.statusText
            })
          )
      })
  })
}

function downloadPDBMeta(pdbidstring) {
  return $.when(
    $.get({
      url:
        "https://www.rcsb.org/pdb/rest/describeMol?structureId=" + pdbidstring,
      dataType: "xml",
      success: function(xml) {
        $(xml)
          .find("structureId")
          .each(function(i, p) {
            var el = $(p)
            var pdbid = el.attr("id")
            //pdb_info[pdbid].macroMolecule = el.find('macroMolecule').attr('name');
            pdb_info[pdbid].chains = el
              .find("chain")
              .map(function(i, d) {
                return d.id
              })
              .toArray()
          })
      }
    }),
    $.get({
      url:
        "https://www.rcsb.org/pdb/rest/describePDB?structureId=" + pdbidstring,
      dataType: "xml",
      success: function(xml) {
        $(xml)
          .find("PDB")
          .each(function(i, p) {
            var pdb = $(p)
            var pdbid = pdb.attr("structureId")
            pdb_info[pdbid].title = pdb.attr("title")
            pdb_info[pdbid].expMethod = pdb.attr("expMethod")
            pdb_info[pdbid].resolution = parseFloat(pdb.attr("resolution"))
            pdb_info[pdbid].deposition_date = pdb.attr("deposition_date")
            pdb_info[pdbid].last_modification_date = pdb.attr(
              "last_modification_date"
            )
            pdb_info[pdbid].structure_authors = pdb.attr("structure_authors")
            pdb_info[pdbid].pubmedId = pdb.attr("pubmedId")
          })
      }
    }),
    $.get({
      url:
        "https://www.rcsb.org/pdb/rest/ligandInfo?structureId=" + pdbidstring,
      dataType: "xml",
      success: function(xml) {
        $(xml)
          .find("structureId")
          .each(function(i, p) {
            var el = $(p)
            var pdbid = el.attr("id")
            var heavy = el
              .find("ligand")
              .toArray()
              .reduce(function(a, b) {
                return +$(a).attr("molecularWeight") >
                  +$(b).attr("molecularWeight")
                  ? a
                  : b
              })
            pdb_info[pdbid].smiles = $(heavy)
              .find("smiles")
              .text()
              .replace(/\\\\/g, "\\")
              .replace("CC(=O)O", "*")
            pdb_info[pdbid].chemId = $(heavy).attr("chemicalID")
            pdb_info[pdbid].chemicalName = $(heavy)
              .find("chemicalName")
              .text()
            pdb_info[pdbid].formula = $(heavy)
              .find("formula")
              .text()
            pdb_info[pdbid].molecularWeight = $(heavy).attr("molecularWeight")
          })
      }
    })
  )
}

async function loadSmiles(pdbid) {
  const { default: SmilesDrawer } = await import("smiles-drawer")

  var canv = $("#smilesCanvas")

  var myDrawer = new SmilesDrawer.Drawer({
    compactDrawing: false,
    height: 300
  })

  canv[0].getContext("2d").clearRect(0, 0, canv.width, canv.height)
  if (canv.data("chemId") !== pdb_info[pdbid].chemId) {
    SmilesDrawer.parse(pdb_info[pdbid].smiles, function(tree) {
      canv.remove()
      canv = $("<canvas>", {
        id: "smilesCanvas"
      }).appendTo($("#smilesDrawing .card-img-top"))
      myDrawer.draw(tree, "smilesCanvas", "light", false)
      canv.data("chemId", pdb_info[pdbid].chemId)
      canv.height("auto")
      canv.width("100%")
    })
  }
}

function loadChemInfo(pdbid) {
  var e = pdb_info[pdbid]
  $("#chem-title").html(e.title)
  $("#chem-date").html(e.deposition_date)
  $("#chem-authors").html(e.structure_authors.split(",")[0] + " et al. ")
  $("#chem-pubmed").attr(
    "href",
    "https://www.ncbi.nlm.nih.gov/pubmed/" + e.pubmedId
  )
  $("#chem-id").html(
    '<a target="_blank" rel="noopener" class="text-secondary" href="https://www.rcsb.org/ligand/' +
      e.chemId +
      '">' +
      e.chemId +
      "</a>"
  )
  $("#chem-form").html(e.formula)
}

//variable to store LiteMol Component Scope which has all the methods
// var liteMolScope
var pdb_info = {}
var getPDBinfo = function(pdbids) {
  let pdbidstring
  if (Array.isArray(pdbids)) {
    for (var i = 0; i < pdbids.length; ++i) {
      pdbids[i] = pdbids[i].toUpperCase()
      pdb_info[pdbids[i]] = {}
    }
    pdbidstring = pdbids.join(",")
  } else {
    pdbids = pdbids.toUpperCase()
    pdb_info[pdbids] = {}
    pdbidstring = pdbids
  }

  downloadPDBMeta(pdbidstring)
    .done(function() {
      var select = $("#pdb_select")
      $(
        pdbids.sort(function(a, b) {
          return pdb_info[a].resolution - pdb_info[b].resolution
        })
      ).each(function(i, d) {
        select.append(
          $("<option>", {
            value: d
          }).html(
            `${d} (${pdb_info[d].resolution} Å, ${
              pdb_info[d].chains.length
            }-chain)`
          )
        )
      })

      initLiteMol("#litemol-viewer", select)
    })
    .fail(function() {
      $(pdbids).each(function(i, d) {
        $("#pdb_select").append(
          $("<option>", {
            value: d
          }).html(d)
        )
      })
      $("#pdb-info").html(
        '<small class="text-muted">failed to retrieve metadata from PDB</small>'
      )
    })
}

// var loadDensity = function() {
//   liteMolScope.LiteMolComponent.loadDensity()
//   $("#density_loader").text("toggle electron density")
//   $("#density_loader").attr("onclick", "ToggleDensity()")
// }

// var ToggleDensity = function() {
//   liteMolScope.LiteMolComponent.toggleDensity()
// }

export default function initPDB(pdbids) {
  getPDBinfo(pdbids)
}

window.initPDB = initPDB

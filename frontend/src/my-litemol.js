import "./css/litemol/LiteMol-plugin-blue.css"
import $ from "jquery"
import LiteMol from "./js/pdb/LiteMol-plugin"

// populated by downloadPDBMeta.success
const pdbInfo = {}

async function loadSmiles(pdbid) {
  let chromophore = pdbInfo[pdbid].chromophore
  if (!chromophore) {
    return
  }
  const _id = pdbInfo[pdbid].chromophore.id
  const url = `https://cdn.rcsb.org/images/ccd/unlabeled/${_id[0]}/${_id}.svg`
  $("#smilesDrawing div").html(
    `<a href="https://www.rcsb.org/ligand/${_id}">
      <img id="smilesImg" src="${url}" alt="Diagram chromophore structure (${_id})">
    </a>`
  )
}

function getPDBbinary(id) {
  return new Promise(function (resolve, reject) {
    $.get(`https://files.rcsb.org/download/${id}.cif`)
      .done((response) => {
        resolve(response)
      })
      .fail((xhr, status, error) => {
        $.get(
          `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}_updated.cif`
        )
          .done((response) => {
            resolve(response)
          })
          .fail((_xhr) =>
            reject(
              new Error({
                status: _xhr.status,
                statusText: _xhr.statusText,
              })
            )
          )
      })
  })
}

function loadChemInfo(pdbid) {
  const e = pdbInfo[pdbid]
  $("#chem-title").html(e.struct.title)
  const d = new Date(e.rcsb_accession_info.deposit_date)
  $("#chem-date").html(
    d.toLocaleDateString("en-US", { year: "numeric", month: "short" })
  )
  $("#chem-authors").html(`${e.audit_author[0].name} et al. `)
  $("#chem-pubmed").attr(
    "href",
    `https://www.ncbi.nlm.nih.gov/pubmed/${e.rcsb_primary_citation.pdbx_database_id_PubMed}`
  )
  if (e.chromophore !== undefined) {
    $("#chem-id").html(
      `<a target="_blank" rel="noopener" class="text-secondary"
      href="https://www.rcsb.org/ligand/${e.chromophore.id}">${e.chromophore.id}</a>`
    )
    $("#chem-form").html(e.chromophore.formula)
  } else {
    $("#chem-id").html("")
  }
}

function initLiteMol(selection, changer) {
  const PluginSpec = LiteMol.Plugin.getDefaultSpecification()
  const { LayoutRegion } = LiteMol.Bootstrap.Components
  const { Components } = LiteMol.Plugin
  PluginSpec.components = [
    Components.Visualization.HighlightInfo(LayoutRegion.Main, true),
    Components.Entity.Current("LiteMol", LiteMol.Plugin.VERSION.number)(
      LayoutRegion.Right,
      true
    ),
    Components.Transform.View(LayoutRegion.Right),
    // Components.Context.Log(LayoutRegion.Bottom, true),
    Components.Context.Overlay(LayoutRegion.Root),
    // Components.Context.Toast(LayoutRegion.Main, true),
    Components.Context.BackgroundTasks(LayoutRegion.Main, true),
  ]

  try {
    const plugin = LiteMol.Plugin.create({
      customSpecification: PluginSpec,
      target: selection,
      viewportBackground: "#fff",
      layoutState: {
        hideControls: true,
        isExpanded: false,
      },
      allowAnalytics: true,
    })

    const dataCache = {}
    changer.change(function () {
      const id = this.value
      if (!Object.prototype.hasOwnProperty.call(dataCache, "id")) {
        dataCache[id] = getPDBbinary(id)
      }
      plugin.clear()
      dataCache[id].then(
        (data) =>
          plugin.loadMolecule({
            data,
            id,
          }),
        (reason) => {
          $(selection).html(
            '<span class="text-danger muted">failed to retrieve molecular structure</span>'
          )
        }
      )
    })

    $("body").on("click", function (e) {
      if ($(".lm-layout-right").length) {
        if ($(e.target).closest("#litemol-viewer").length === 0) {
          plugin.setLayoutState({
            hideControls: true,
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
    .change(function () {
      // var id = this.value
      $("#pdb-out-link").attr(
        "href",
        `https://www.rcsb.org/structure/${this.value}`
      )
      if (pdbInfo[this.value]) {
        loadSmiles(this.value)
        loadChemInfo(this.value)
      }
    })
    .trigger("change")
}

function downloadPDBMeta(pdbIds) {
  return $.when(
    $.post({
      url: "https://data.rcsb.org/graphql",
      contentType: "application/json",
      data: JSON.stringify({
        query: `{
          entries(entry_ids:${pdbIds}) {
            entry {
              id
            }
            struct {
              title
            }
            rcsb_entry_info {
              experimental_method
              resolution_combined
            }
            rcsb_accession_info {
              deposit_date
              revision_date
            }
            audit_author {
              name
            }
            rcsb_primary_citation {
              pdbx_database_id_PubMed
            }
            polymer_entities {
              chem_comp_nstd_monomers {
                chem_comp {
                  id
                  type
                  formula_weight
                  name
                  formula
                }
              }
            }
            nonpolymer_entities {
              nonpolymer_comp {
                chem_comp {
                  id
                  type
                  formula_weight
                  name
                  formula
                }
              }
            }
          }
        }`,
      }),
      success: function ({ data }) {
        data.entries.forEach((entry) => {
          pdbInfo[entry.entry.id] = entry
          let chromo = null

          if (entry.polymer_entities && entry.polymer_entities.length > 0) {
            chromo = entry.polymer_entities[0].chem_comp_nstd_monomers
          } else if (
            entry.nonpolymer_entities &&
            entry.nonpolymer_entities.length > 0
          ) {
            chromo = entry.nonpolymer_entities[0].nonpolymer_comp
          }

          if (Array.isArray(chromo)) {
            chromo = chromo.reduce((a, b) =>
              a.chem_comp.formula_weight > b.chem_comp.formula_weight ? a : b
            )
          }

          if (chromo && chromo.chem_comp) {
            pdbInfo[entry.entry.id].chromophore = { ...chromo.chem_comp }
          }

          if (
            entry.rcsb_entry_info &&
            entry.rcsb_entry_info.resolution_combined &&
            entry.rcsb_entry_info.resolution_combined.length > 0
          ) {
            pdbInfo[entry.entry.id].resolution =
              entry.rcsb_entry_info.resolution_combined[0]
          }
        })
      },
    })
  )
}

const getPDBinfo = function (pdbIds) {
  downloadPDBMeta(`["${pdbIds.join('","')}"]`)
    .done(() => {
      const select = $("#pdb_select")
      pdbIds.sort((a, b) =>
        pdbInfo[a].resolution > pdbInfo[b].resolution ? 1 : -1
      )
      pdbIds.forEach((id) => {
        select.append(
          $("<option>", {
            value: id,
          }).html(`${id} (${pdbInfo[id].resolution} Å)`)
        )
      })

      initLiteMol("#litemol-viewer", select)
    })
    .fail(function () {
      let html =
        '<div><p><small class="text-muted">Failed to retrieve metadata from PDB!</small></p>'
      html += "Please look for these IDs at RSCB PDB: &nbsp;"
      pdbIds.forEach((_id) => {
        html += `<a href="https://www.rcsb.org/structure/${_id}">${_id}</a>, `
      })
      html += "</div>"
      $("#protein-structure").html(html).removeClass("row")
    })
}

export default function initPDB(pdbids) {
  getPDBinfo(pdbids)
}

window.initPDB = initPDB

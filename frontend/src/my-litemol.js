// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/jquery-ajax-sentry.js" // Track jQuery AJAX errors

import "./css/litemol/LiteMol-plugin-blue.css"

// Import UMD bundle - it sets window.LiteMol global
import "./js/pdb/LiteMol-plugin"
const LiteMol = window.LiteMol

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE.currentBundle = "litemol"

// populated by downloadPDBMeta.success
const pdbInfo = {}

/**
 * Load and display chromophore structure diagram from RCSB.
 *
 * @param {string} pdbid - PDB identifier
 */
function loadSmiles(pdbid) {
  const chromophore = pdbInfo[pdbid]?.chromophore
  if (!chromophore) return

  // Escape HTML to prevent XSS
  const chromoId = $("<div>").text(chromophore.id).html()
  const svgUrl = `https://cdn.rcsb.org/images/ccd/unlabeled/${chromoId[0]}/${chromoId}.svg`

  $("#smilesDrawing div").html(
    `<a href="https://www.rcsb.org/ligand/${chromoId}">
      <img id="smilesImg" src="${svgUrl}" alt="Chromophore structure diagram (${chromoId})">
    </a>`
  )
}

/**
 * Fetch PDB structure file from multiple mirror sources with cascading fallback.
 *
 * Tries sources in order of preference:
 * 1. wwPDB (official worldwide PDB, recommended by RCSB docs)
 * 2. RCSB (US mirror, byte-identical to wwPDB)
 * 3. EBI/PDBe (European mirror, compatible format)
 *
 * @param {string} id - PDB identifier (e.g., '6GP0')
 * @returns {Promise<string>} CIF format structure data
 * @throws {Error} If all endpoints fail
 * @see https://www.rcsb.org/docs/programmatic-access/file-download-services
 */
async function getPDBbinary(id) {
  const TIMEOUT = 15000 // 15 second timeout per request

  const endpoints = [
    {
      name: "wwPDB",
      url: `https://files.wwpdb.org/download/${id}.cif`,
    },
    {
      name: "RCSB",
      url: `https://files.rcsb.org/download/${id}.cif`,
    },
    {
      name: "EBI",
      url: `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}_updated.cif`,
    },
  ]

  const errors = []

  // Try each endpoint in sequence
  for (let i = 0; i < endpoints.length; i++) {
    const { name, url } = endpoints[i]
    const isLastEndpoint = i === endpoints.length - 1

    try {
      const response = await $.ajax({ url, timeout: TIMEOUT })
      return response
    } catch (error) {
      // Track the error
      errors.push({ endpoint: name, error })

      // Log to Sentry for monitoring
      if (window.Sentry) {
        window.Sentry.addBreadcrumb({
          message: `PDB fetch failed: ${name}`,
          data: {
            pdbId: id,
            endpoint: name,
            error: error.statusText,
            attemptNumber: i + 1,
            totalEndpoints: endpoints.length,
          },
          level: "warning",
        })
      }

      // If this was the last endpoint, throw with all error details
      if (isLastEndpoint) {
        const errorSummary = errors
          .map((e) => `${e.endpoint}: ${e.error.statusText || "Network error"}`)
          .join("; ")

        throw new Error(`Failed to fetch PDB ${id} from all sources. Errors: ${errorSummary}`)
      }
    }
  }
}

/**
 * Load and display chemical information for a PDB entry.
 *
 * @param {string} pdbid - PDB identifier
 */
function loadChemInfo(pdbid) {
  const entry = pdbInfo[pdbid]
  if (!entry) return

  // Use .text() for safe insertion of untrusted content
  $("#chem-title").text(entry.struct?.title || "Unknown")

  const depositDate = entry.rcsb_accession_info?.deposit_date
  if (depositDate) {
    const date = new Date(depositDate)
    $("#chem-date").html(
      date.toLocaleDateString("en-US", {
        year: "numeric",
        month: "short",
      })
    )
  }

  // Safe access to potentially missing author data
  const firstAuthor = entry.audit_author?.[0]?.name
  if (firstAuthor) {
    $("#chem-authors").text(`${firstAuthor} et al. `)
  }

  // Safe access to PubMed ID
  const pubmedId = entry.rcsb_primary_citation?.pdbx_database_id_PubMed
  if (pubmedId) {
    $("#chem-pubmed").attr("href", `https://www.ncbi.nlm.nih.gov/pubmed/${pubmedId}`)
  } else {
    $("#chem-pubmed").removeAttr("href")
  }

  if (entry.chromophore) {
    // Escape HTML to prevent XSS
    const chromoId = $("<div>").text(entry.chromophore.id).html()
    $("#chem-id").html(
      `<a target="_blank" rel="noopener" class="text-secondary"
      href="https://www.rcsb.org/ligand/${chromoId}">${chromoId}</a>`
    )
    $("#chem-form").text(entry.chromophore.formula)
  } else {
    $("#chem-id").html("")
  }
}

/**
 * Initialize LiteMol molecular visualization plugin.
 *
 * @param {string} selection - CSS selector for the container element
 * @param {jQuery} changer - jQuery object for the PDB selector dropdown
 */
function initLiteMol(selection, changer) {
  const PluginSpec = LiteMol.Plugin.getDefaultSpecification()
  const { LayoutRegion } = LiteMol.Bootstrap.Components
  const { Components } = LiteMol.Plugin
  PluginSpec.components = [
    Components.Visualization.HighlightInfo(LayoutRegion.Main, true),
    Components.Entity.Current("LiteMol", LiteMol.Plugin.VERSION.number)(LayoutRegion.Right, true),
    Components.Transform.View(LayoutRegion.Right),
    Components.Context.Overlay(LayoutRegion.Root),
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

    // Cache PDB data promises to avoid redundant fetches
    const dataCache = new Map()
    let currentRequest = null

    changer.change(async function () {
      const id = this.value
      const requestId = Symbol("request")
      currentRequest = requestId

      // Get or create cached promise
      if (!dataCache.has(id)) {
        dataCache.set(id, getPDBbinary(id))
      }

      plugin.clear()

      try {
        const data = await dataCache.get(id)

        // Only load if this is still the active request
        if (currentRequest === requestId) {
          plugin.loadMolecule({ data, id })
        }
      } catch (error) {
        // Only show error if this is still the active request
        if (currentRequest === requestId) {
          $(selection).html(
            '<span class="text-danger muted">Failed to retrieve molecular structure. Please refresh.</span>'
          )
          if (window.Sentry) {
            window.Sentry.captureException(error, {
              tags: { pdbId: id, component: "litemol" },
            })
          }
        }
      }
    })

    // Close side panel when clicking outside (use namespace to prevent leaks)
    $("body")
      .off("click.litemol")
      .on("click.litemol", (e) => {
        if ($(".lm-layout-right").length && $(e.target).closest("#litemol-viewer").length === 0) {
          plugin.setLayoutState({ hideControls: true })
        }
      })
  } catch (err) {
    if (window.Sentry) {
      window.Sentry.captureException(err, {
        tags: { component: "litemol-init" },
      })
    }
  }

  // Update external link and load metadata when selection changes
  changer
    .change(function () {
      const id = this.value
      $("#pdb-out-link").attr("href", `https://www.rcsb.org/structure/${id}`)

      if (pdbInfo[id]) {
        loadSmiles(id)
        loadChemInfo(id)
      }
    })
    .trigger("change")
}

/**
 * Download PDB metadata from RCSB GraphQL API.
 *
 * @param {string} pdbIds - JSON array string of PDB IDs (e.g., '["6GP0","1GFL"]')
 * @returns {Promise} jQuery promise that resolves when metadata is loaded
 */
function downloadPDBMeta(pdbIds) {
  // Check cache first (7-day TTL per RCSB recommendations)
  const cacheKey = `pdb_meta_${pdbIds}`
  const CACHE_TTL = 7 * 24 * 60 * 60 * 1000 // 7 days in milliseconds

  try {
    const cached = localStorage.getItem(cacheKey)
    if (cached) {
      const { timestamp, data } = JSON.parse(cached)
      if (Date.now() - timestamp < CACHE_TTL) {
        // Restore cached data to pdbInfo
        Object.assign(pdbInfo, data)
        return Promise.resolve()
      }
    }
  } catch (error) {
    // Ignore cache errors, fall through to fetch
    console.warn("Cache read failed:", error)
  }

  return $.post({
    url: "https://data.rcsb.org/graphql",
    contentType: "application/json",
    data: JSON.stringify({
      query: `{
        entries(entry_ids:${pdbIds}) {
          entry { id }
          struct { title }
          rcsb_entry_info {
            experimental_method
            resolution_combined
          }
          rcsb_accession_info {
            deposit_date
            revision_date
          }
          audit_author { name }
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
  }).then((response) => {
    // GraphQL always returns 200 OK, check for errors in response body
    if (response.errors) {
      const errorMsg = response.errors.map((e) => e.message).join("; ")
      throw new Error(`GraphQL errors: ${errorMsg}`)
    }

    if (!response.data?.entries) {
      throw new Error("No data returned from GraphQL API")
    }
    const { data } = response
    const fetchedData = {}

    data.entries.forEach((entry) => {
      const entryId = entry.entry.id
      pdbInfo[entryId] = entry
      fetchedData[entryId] = entry

      // Extract chromophore (largest component by molecular weight)
      let chromo = null
      if (entry.polymer_entities?.length > 0) {
        chromo = entry.polymer_entities[0].chem_comp_nstd_monomers
      } else if (entry.nonpolymer_entities?.length > 0) {
        chromo = entry.nonpolymer_entities[0].nonpolymer_comp
      }

      if (Array.isArray(chromo) && chromo.length > 0) {
        chromo = chromo.reduce((a, b) =>
          a.chem_comp.formula_weight > b.chem_comp.formula_weight ? a : b
        )
      }

      if (chromo?.chem_comp) {
        pdbInfo[entryId].chromophore = { ...chromo.chem_comp }
      }

      // Extract resolution if available
      const resolutions = entry.rcsb_entry_info?.resolution_combined
      if (resolutions?.length > 0) {
        pdbInfo[entryId].resolution = resolutions[0]
      }
    })

    // Cache the successful result
    try {
      localStorage.setItem(
        cacheKey,
        JSON.stringify({
          timestamp: Date.now(),
          data: fetchedData,
        })
      )
    } catch (error) {
      // Ignore cache write errors (quota exceeded, etc.)
      console.warn("Cache write failed:", error)
    }
  })
}

/**
 * Fetch PDB metadata and initialize LiteMol viewer.
 *
 * @param {string[]} pdbIds - Array of PDB identifiers
 */
async function getPDBinfo(pdbIds) {
  try {
    await downloadPDBMeta(`["${pdbIds.join('","')}"]`)

    const select = $("#pdb_select")

    // Sort by resolution (best first)
    pdbIds.sort((a, b) => (pdbInfo[a].resolution > pdbInfo[b].resolution ? 1 : -1))

    // Populate dropdown
    pdbIds.forEach((id) => {
      select.append($("<option>", { value: id }).html(`${id} (${pdbInfo[id].resolution} Å)`))
    })

    initLiteMol("#litemol-viewer", select)
  } catch (error) {
    // Log error to Sentry for monitoring
    if (window.Sentry) {
      window.Sentry.captureException(error, {
        tags: { component: "pdb-metadata", pdbIds: pdbIds.join(",") },
      })
    }
    console.error("Failed to load PDB metadata:", error)

    const links = pdbIds
      .map((id) => `<a href="https://www.rcsb.org/structure/${id}">${id}</a>`)
      .join(", ")

    $("#protein-structure")
      .html(
        `<div>
          <p class="text-danger muted">Failed to retrieve metadata from PDB!</p>
          <p>You may view these PDB IDs directly at RCSB: ${links}</p>
        </div>`
      )
      .removeClass("row")
  }
}

export default function initPDB(pdbids) {
  getPDBinfo(pdbids)
}

window.initPDB = initPDB

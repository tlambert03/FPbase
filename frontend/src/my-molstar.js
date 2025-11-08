// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/jquery-ajax-sentry.js" // Track jQuery AJAX errors

// Import Mol* (pdbe-molstar) plugin
import { PDBeMolstarPlugin } from "pdbe-molstar/lib/viewer"
import "pdbe-molstar/build/pdbe-molstar.css"

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE.currentBundle = "molstar"

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
 * Get PDB structure URL with cascading fallback.
 *
 * Returns the first available URL from multiple mirror sources.
 * Tries in order: wwPDB → RCSB → EBI
 *
 * @param {string} id - PDB identifier (e.g., '6GP0')
 * @returns {Promise<string>} URL to fetch CIF structure data
 * @throws {Error} If all endpoints fail
 * @see https://www.rcsb.org/docs/programmatic-access/file-download-services
 */
async function getPDBUrl(id) {
  const TIMEOUT = 5000 // 5 second timeout per request

  // Try URLs in order of preference
  const urls = [
    `https://files.wwpdb.org/download/${id}.cif`,
    `https://files.rcsb.org/download/${id}.cif`,
    `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}_updated.cif`,
  ]

  // Test each URL with a HEAD request to see if it's available
  for (const url of urls) {
    try {
      await $.ajax({ url, method: "HEAD", timeout: TIMEOUT })
      return url // This URL works!
    } catch (error) {}
  }

  // All URLs failed - return the first one and let pdbe-molstar try
  // (it will show its own error message)
  return urls[0]
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
 * Initialize Mol* (pdbe-molstar) molecular visualization plugin.
 *
 * @param {string} selection - CSS selector for the container element
 * @param {jQuery} changer - jQuery object for the PDB selector dropdown
 */
function initMolstar(selection, changer) {
  const viewerContainer = document.querySelector(selection)
  if (!viewerContainer) {
    console.error(`Mol* viewer container not found: ${selection}`)
    return
  }

  // Create plugin instance
  const plugin = new PDBeMolstarPlugin()

  let currentPluginInstance = null

  changer.change(async function () {
    const id = this.value

    // Clear previous structure
    if (currentPluginInstance) {
      // Clear the container for re-rendering
      viewerContainer.innerHTML = ""
    }

    try {
      // Get a working URL with cascading fallback
      const url = await getPDBUrl(id)

      // Render the plugin with custom data URL
      const options = {
        customData: { url: url, format: "cif" },
        bgColor: { r: 255, g: 255, b: 255 },
        hideControls: true,
        sequencePanel: true,
        hideStructure: ["het", "water", "carbs"],
      }

      await plugin.render(viewerContainer, options)

      // Hide the axis helper (XYZ indicator) at bottom left
      plugin.plugin.canvas3d?.setProps({
        camera: {
          helper: { axes: { name: "off", params: {} } },
        },
      })

      // Monitor layout state changes to enforce controls visibility rules
      // In standard (non-expanded) mode, controls must remain hidden
      plugin.plugin.layout.events.updated.subscribe(() => {
        const state = plugin.plugin.layout.state
        if (!state.isExpanded && state.showControls) {
          // User exited expanded mode but controls are visible - hide them
          plugin.canvas.toggleControls(false)
        }
      })

      currentPluginInstance = plugin
    } catch (error) {
      $(selection).html(
        '<span class="text-danger muted">Failed to retrieve molecular structure. Please refresh.</span>'
      )
      if (window.Sentry) {
        window.Sentry.captureException(error, {
          tags: { pdbId: id, component: "molstar" },
        })
      }
    }
  })

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
 * Fetch PDB metadata and initialize Mol* viewer.
 *
 * @param {string[]} pdbIds - Array of PDB identifiers
 */
async function getPDBinfo(pdbIds) {
  try {
    await downloadPDBMeta(`["${pdbIds.join('","')}"]`)

    const select = $("#pdb_select")

    // Sort by resolution (best first), handling missing resolution data
    // Entries without resolution (e.g., NMR structures) are sorted to the end
    pdbIds.sort((a, b) => {
      const resA = pdbInfo[a]?.resolution
      const resB = pdbInfo[b]?.resolution

      // If both have resolution, sort by value (lower is better)
      if (resA !== undefined && resB !== undefined) {
        return resA - resB
      }
      // If only A has resolution, it comes first
      if (resA !== undefined) return -1
      // If only B has resolution, it comes first
      if (resB !== undefined) return 1
      // If neither has resolution, maintain original order
      return 0
    })

    // Populate dropdown
    pdbIds.forEach((id) => {
      const resolution = pdbInfo[id]?.resolution
      const displayText =
        resolution !== undefined ? `${id} (${resolution.toFixed(2)} Å)` : `${id} (resolution N/A)`
      select.append($("<option>", { value: id }).html(displayText))
    })

    initMolstar("#molstar-viewer", select)
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

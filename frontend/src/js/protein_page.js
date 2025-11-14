const $ = window.jQuery // jQuery loaded from CDN

window.initSnapGene = (protein, selection) => {
  // Search SnapGene's plasmid database for this protein
  $.get({
    url: "https://www.snapgene.com/api/plasmids/search",
    data: { string: protein },
    crossDomain: true,
    success: (data) => {
      if (!data || !Array.isArray(data)) return

      // Filter to only the Fluorescent Protein Genes & Plasmids set
      const fpSet = data.find(
        (set) => set.setData?.set === "fluorescent_protein_genes_and_plasmids"
      )
      if (!fpSet || !fpSet.sequences) return

      // Filter plasmids to only include relevant matches:
      // - Exact match (e.g., "EGFP")
      // - Plasmid vectors (e.g., "pEGFP", "pEGFP-N1", "pEGFP-C2")
      const isPlasmidVector = new RegExp(`^p${protein}(-[CN]?\\d+)?$`, "i")

      const plasmids = fpSet.sequences
        .filter((plasmid) => {
          const name = plasmid.plasmidName
          return name === protein || isPlasmidVector.test(name)
        })
        .map((plasmid) => ({
          name: plasmid.plasmidName,
          url: `https://www.snapgene.com/plasmids/${plasmid.setID}/${plasmid.plasmidID}`,
        }))

      if (plasmids.length > 0) {
        const label = plasmids.length === 1 ? "SnapGene plasmid: " : "SnapGene plasmids: "
        const $li = $("<li>").text(label).appendTo($(selection))

        plasmids.forEach((plasmid, index) => {
          $li.append(
            $("<a>", {
              href: plasmid.url,
              target: "_blank",
              rel: "noopener",
            }).text(plasmid.name)
          )
          if (index !== plasmids.length - 1) {
            $li.append(", ")
          }
        })
      }
    },
    error: () => {
      // Silently fail - this is a nice-to-have feature
      console.debug("SnapGene search unavailable")
    },
  })
}

// Toggle between relative and root mutations display
window.toggleMutationDisplay = (event) => {
  event.preventDefault()
  const clickedElement = event.target
  const isRelParent = clickedElement.classList.contains("mut-rel-parent")

  const classToHide = isRelParent ? "mut-rel-parent" : "mut-rel-root"
  const classToShow = isRelParent ? "mut-rel-root" : "mut-rel-parent"

  document.querySelectorAll(`.${classToHide}`).forEach((el) => {
    el.classList.add("hidden")
  })

  document.querySelectorAll(`.${classToShow}`).forEach((el) => {
    el.classList.remove("hidden")
  })
}

// Initialize lineage chart with IntersectionObserver for lazy loading
window.initLineageChart = (slug) => {
  const element = document.querySelector(".lineage")
  if (!element) return

  // Use IntersectionObserver to only load when scrolled into view
  const observer = new IntersectionObserver(
    (entries) => {
      entries.forEach((entry) => {
        if (entry.isIntersecting) {
          loadLineageChart(slug, element)
          observer.unobserve(entry.target)
        }
      })
    },
    {
      rootMargin: "200px", // Load when 200px from entering viewport
    }
  )

  observer.observe(element)
}

function loadLineageChart(slug, element) {
  // Wait for FPBASE to be ready, then load D3 charts
  function initLineage() {
    window.FPBASE.loadD3Charts()
      .then(({ LineageChart }) => {
        $.getJSON(`/ajax/lineage/${slug}/`, (data) => {
          window._lastLineageData = data
          if (!data.children || data.children.length === 0) {
            return
          }
          const linchart = LineageChart({ slug }).data(data)
          const lineage = d3.select(element)
          lineage.call(linchart)
          linchart.duration(200)
          if (slug && $(".lineage-wrapper").length) {
            const node = d3.select(`#node_${slug}`)
            if (!node.empty()) {
              const slugpos = node.datum().y
              $(".lineage-wrapper").scrollLeft(slugpos - (window.innerWidth - 30) / 3)
            }
          }
        })
      })
      .catch((error) => {
        console.error("Error loading D3 charts:", error)
        $(element).html(
          '<div class="alert alert-danger">Failed to load lineage viewer. Please refresh the page.</div>'
        )
      })
  }

  if (window.FPBASE && window.FPBASE.loadD3Charts) {
    initLineage()
  } else {
    window.addEventListener("fpbase:ready", initLineage)
  }
}

import $ from 'jquery'

window.initSnapGene = (protein, selection) => {
  // Search SnapGene's plasmid database for this protein
  $.get({
    url: 'https://www.snapgene.com/api/plasmids/search',
    data: { string: protein },
    crossDomain: true,
    success: (data) => {
      if (!data || !Array.isArray(data)) return

      // Filter to only the Fluorescent Protein Genes & Plasmids set
      const fpSet = data.find(
        (set) => set.setData?.set === 'fluorescent_protein_genes_and_plasmids'
      )
      if (!fpSet || !fpSet.sequences) return

      // Filter plasmids to only include relevant matches:
      // - Exact match (e.g., "EGFP")
      // - Plasmid vectors (e.g., "pEGFP", "pEGFP-N1", "pEGFP-C2")
      const isPlasmidVector = new RegExp(`^p${protein}(-[CN]?\\d+)?$`, 'i')

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
        const label = plasmids.length === 1 ? 'SnapGene plasmid: ' : 'SnapGene plasmids: '
        const $li = $('<li>').text(label).appendTo($(selection))

        plasmids.forEach((plasmid, index) => {
          $li.append(
            $('<a>', {
              href: plasmid.url,
              target: '_blank',
              rel: 'noopener',
            }).text(plasmid.name)
          )
          if (index !== plasmids.length - 1) {
            $li.append(', ')
          }
        })
      }
    },
    error: () => {
      // Silently fail - this is a nice-to-have feature
      console.debug('SnapGene search unavailable')
    },
  })
}

import $ from "jquery"

window.initSnapGene = function(protein, selection) {
  const snaptemplate =
    "https://www.snapgene.com/resources/plasmid_files/fluorescent_protein_genes_and_plasmids/*/"
  $.get({
    url: "https://fpbase.s3.amazonaws.com/snapgene.xml",
    crossDomain: true,
    success: function(data) {
      /* handle data here */
      const hits = []
      $(data)
        .find(`[name*="${protein}"]`)
        .each(function() {
          let name = $(this).attr("name")
          while (name.charAt(0) === "p") {
            name = name.substr(1)
          }
          $.each(
            ["-N2", "-N1", "-N", "-C2", "-C1", "-C", "-1", "-B", "-S1", "-PRL"],
            function(index, val) {
              name = name.replace(val, "")
            }
          )
          if (name === protein || name === `${protein}1`) {
            hits.push([$(this).attr("name"), $(this).attr("url")])
          }
        })

      if (hits.length) {
        let t = "SnapGene links: "
        if (hits.length === 1) {
          t = "SnapGene link: "
        }
        const $li = $("<li>")
          .text(t)
          .appendTo($(selection))
        $.each(hits, function(index, val) {
          $li.append(
            $("<a>", { href: snaptemplate.replace("*", val[1]) }).text(val[0])
          )
          if (index !== hits.length - 1) {
            $li.append(", ")
          }
        })
      }
    },
  })
}

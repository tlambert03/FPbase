import { icon } from "./icons.js"

const $ = window.jQuery // jQuery loaded from CDN
import "./detect-touch" // adds window.USER_IS_TOUCHING = true; after touch event.

// Helper to wait for Bootstrap plugins to be available
function waitForBootstrap(callback) {
  if (typeof $.fn.tooltip !== "undefined" && typeof $.fn.popover !== "undefined") {
    callback()
  } else {
    // Retry after a short delay
    setTimeout(() => waitForBootstrap(callback), 50)
  }
}

window.mobilecheck = () => {
  var check = false
  ;((a) => {
    if (
      /(android|bb\d+|meego).+mobile|avantgo|bada\/|blackberry|blazer|compal|elaine|fennec|hiptop|iemobile|ip(hone|od)|iris|kindle|lge |maemo|midp|mmp|mobile.+firefox|netfront|opera m(ob|in)i|palm( os)?|phone|p(ixi|re)\/|plucker|pocket|psp|series(4|6)0|symbian|treo|up\.(browser|link)|vodafone|wap|windows ce|xda|xiino/i.test(
        a
      ) ||
      /1207|6310|6590|3gso|4thp|50[1-6]i|770s|802s|a wa|abac|ac(er|oo|s-)|ai(ko|rn)|al(av|ca|co)|amoi|an(ex|ny|yw)|aptu|ar(ch|go)|as(te|us)|attw|au(di|-m|r |s )|avan|be(ck|ll|nq)|bi(lb|rd)|bl(ac|az)|br(e|v)w|bumb|bw-(n|u)|c55\/|capi|ccwa|cdm-|cell|chtm|cldc|cmd-|co(mp|nd)|craw|da(it|ll|ng)|dbte|dc-s|devi|dica|dmob|do(c|p)o|ds(12|-d)|el(49|ai)|em(l2|ul)|er(ic|k0)|esl8|ez([4-7]0|os|wa|ze)|fetc|fly(-|_)|g1 u|g560|gene|gf-5|g-mo|go(\.w|od)|gr(ad|un)|haie|hcit|hd-(m|p|t)|hei-|hi(pt|ta)|hp( i|ip)|hs-c|ht(c(-| |_|a|g|p|s|t)|tp)|hu(aw|tc)|i-(20|go|ma)|i230|iac( |-|\/)|ibro|idea|ig01|ikom|im1k|inno|ipaq|iris|ja(t|v)a|jbro|jemu|jigs|kddi|keji|kgt( |\/)|klon|kpt |kwc-|kyo(c|k)|le(no|xi)|lg( g|\/(k|l|u)|50|54|-[a-w])|libw|lynx|m1-w|m3ga|m50\/|ma(te|ui|xo)|mc(01|21|ca)|m-cr|me(rc|ri)|mi(o8|oa|ts)|mmef|mo(01|02|bi|de|do|t(-| |o|v)|zz)|mt(50|p1|v )|mwbp|mywa|n10[0-2]|n20[2-3]|n30(0|2)|n50(0|2|5)|n7(0(0|1)|10)|ne((c|m)-|on|tf|wf|wg|wt)|nok(6|i)|nzph|o2im|op(ti|wv)|oran|owg1|p800|pan(a|d|t)|pdxg|pg(13|-([1-8]|c))|phil|pire|pl(ay|uc)|pn-2|po(ck|rt|se)|prox|psio|pt-g|qa-a|qc(07|12|21|32|60|-[2-7]|i-)|qtek|r380|r600|raks|rim9|ro(ve|zo)|s55\/|sa(ge|ma|mm|ms|ny|va)|sc(01|h-|oo|p-)|sdk\/|se(c(-|0|1)|47|mc|nd|ri)|sgh-|shar|sie(-|m)|sk-0|sl(45|id)|sm(al|ar|b3|it|t5)|so(ft|ny)|sp(01|h-|v-|v )|sy(01|mb)|t2(18|50)|t6(00|10|18)|ta(gt|lk)|tcl-|tdg-|tel(i|m)|tim-|t-mo|to(pl|sh)|ts(70|m-|m3|m5)|tx-9|up(\.b|g1|si)|utst|v400|v750|veri|vi(rg|te)|vk(40|5[0-3]|-v)|vm40|voda|vulc|vx(52|53|60|61|70|80|81|83|85|98)|w3c(-| )|webc|whit|wi(g |nc|nw)|wmlb|wonu|x700|yas-|your|zeto|zte-/i.test(
        a.substr(0, 4)
      )
    )
      check = true
  })(navigator.userAgent || navigator.vendor || window.opera)
  return check
}

if (document.getElementById("comparison-slider")) {
  if (window.USER_IS_TOUCHING) {
    $("#comparison-toggle").click(() => {
      $("#comparison-slider").toggleClass("hover-effect")
    })
    $(document).on("click", (e) => {
      if (!document.getElementById("comparison-slider").contains(e.target)) {
        $("#comparison-slider").removeClass("hover-effect")
      }
    })
  } else {
    $("#comparison-slider").hover(
      () => {
        $("#comparison-slider").addClass("hover-effect")
      },
      () => {
        $("#comparison-slider").removeClass("hover-effect")
      }
    )
  }
}

$(".custom-file-input").on("change", function () {
  var fileName = $(this).val().split("\\").pop()
  if (fileName === "") {
    fileName = "Choose file"
  }
  $(this).next(".custom-file-label").addClass("selected").html(fileName)
})

$(() => {
  var $quote = $(".protein .name:first")
  var $numChar = $quote.text().length

  if ($numChar >= 40) {
    $quote.css("font-size", "28px")
  } else if ($numChar >= 35) {
    $quote.css("font-size", "32px")
  } else if ($numChar >= 28) {
    $quote.css("font-size", "36px")
  } else {
  }
})

$(".form-group").removeClass("row")

//Navbar Scroll Event.
// comment this out to remove "navbar hiding" when the user scrolls down
// var lastScrollTop = 0;
// var navbar        = $('.navbar');
// $(window).scroll(function(event){
//    var st = $(this).scrollTop();
//    if (st > lastScrollTop){
//        navbar.addClass('navbar-scroll-custom');
//    } else {
//       navbar.removeClass('navbar-scroll-custom');
//    }
//    lastScrollTop = st;
// });

/////////////// COMPARISON SLIDER ///////////

function populate_comparison_tab(comparison_set) {
  var $ul = $("#comparison-slider ul.comparison-list")
  //$ul.empty();
  if (comparison_set.length) {
    //var token = $("#csrfform input").val()
    const currents = $(".comparison-list li")
      .map((_i, v) => $(v).attr("value"))
      .toArray()
    $.each(comparison_set, (_index, val) => {
      if (currents.indexOf(val.slug) >= 0) {
        return true
      }
      let exemstring
      if (val.exMax && val.emMax) {
        exemstring =
          "Ex/Em &lambda;: &nbsp;<strong>" +
          (val.exMax || "") +
          "</strong> / <strong>" +
          val.emMax +
          "</strong>"
      }
      let ecqystring
      if (val.ec && val.qy) {
        const ec = val.ec.toLocaleString()
        ecqystring =
          "<br>EC: <strong>" +
          ec +
          "</strong>&nbsp;&nbsp;&nbsp;QY: <strong>" +
          (val.qy || "") +
          "</strong>"
      }
      var widget = $("<li>", { class: "comparison-item", value: val.slug })
        .append(
          $("<a>", {
            href: `/protein/${val.slug}`,
            style: `color: ${val.color}`,
          }).html(val.name)
        )
        .append($("<p>").html((exemstring || "") + (ecqystring || "")))

      // Only append spectrum image if the protein has spectra
      if (val.spectra && JSON.parse(val.spectra).length > 0) {
        widget.append(
          $("<img>", {
            src: `/spectra_img/${val.slug}.svg?xlim=400,700&fill=1&xlabels=0`,
            class: "img-fluid spectrum-svg",
            onerror: "this.style.display='none'",
            alt: `${val.name} spectrum`,
          })
        )
      }

      widget.append(
        $("<button>", {
          class: "comparison-btn remove-protein",
          "data-op": "remove",
          "data-object": val.slug,
          "data-action-url": "/ajax/comparison/",
        }).html("&times;")
      )
      widget.appendTo($ul)
    })
    $("#clearbutton").show()
    if (comparison_set.length === 1) {
      $("#compare-link a").hide()
      $("#compare-link div.msg").text("Add at least two proteins...")
    } else {
      $("#compare-link a").show()
      $("#compare-link div.msg").text("")
    }
  } else {
    $("#clearbutton").hide()
    $("#compare-link a").hide()
    $("#compare-link div.msg").text("Nothing added to comparison...")
  }
}

function handle_comparison_button(e) {
  var button = $(this)
  e.preventDefault()
  $.ajax({
    // create an AJAX call...
    data: {
      object: button.data("object"),
      csrfmiddlewaretoken: window.CSRF_TOKEN,
      operation: button.data("op"),
    },
    type: "POST",
    url: button.attr("data-action-url"),
    dataType: "json",
    success: (response) => {
      // on success..
      populate_comparison_tab(response.comparison_set)
    },
  })
  if ($(this).data("op") === "remove") {
    $(this).closest("li").remove()
  } else if ($(this).data("op") === "clear") {
    $(".comparison-list").empty()
  }
  if ($(this).data("flash")) {
    $("#comparison-toggle").fadeTo(30, 0.3, function () {
      $(this).fadeTo(200, 1.0)
    })
  }
  return false
}

$(document).on("click", ".comparison-btn", handle_comparison_button)
$(".comparison-btn").on("click", handle_comparison_button)

$(() => {
  if (document.getElementById("comparison-slider")) {
    $.getJSON("/ajax/comparison/").then((d) => {
      populate_comparison_tab(d.comparison_set)
    })
  }
})

/////////////////. Spectra Image URL Builder

// Helper to wait for select2 plugin to be available
function waitForSelect2(callback) {
  if (typeof $.fn.select2 !== "undefined") {
    callback()
  } else {
    // Retry after a short delay
    setTimeout(() => waitForSelect2(callback), 50)
  }
}

// Only initialize Select2 if the element exists on the page and not already initialized
if (
  document.getElementById("proteinSlug") &&
  !$("#proteinSlug").hasClass("select2-hidden-accessible")
) {
  waitForSelect2(() => {
    $("#proteinSlug").select2({
      theme: "bootstrap-5",
      width: "80%",
      ajax: {
        theme: "bootstrap-5",
        containerCssClass: ":all:",
        width: "auto",

        url: "/autocomplete-protein",
        dataType: "json",
        cache: true,
        data: (params) => {
          var query = {
            q: params.term,
            type: "spectra",
            page: params.page,
            _type: params._type,
          }
          return query
        },
        processResults: (data) => {
          // Ensure all results have an id property (fixes FPBASE-5H3)
          // The autocomplete endpoint returns results, but we need to ensure
          // each item has an id that can be safely converted to string
          if (data.results) {
            data.results = data.results.map((item) => {
              // If item doesn't have an id, use text as fallback
              if (!item.id && item.text) {
                item.id = item.text
              }
              return item
            })
          }
          return data
        },
      },
      // Safely handle cases where id might still be undefined
      templateResult: (item) => {
        if (!item.id) {
          return item.text
        }
        return item.text || item.id
      },
      templateSelection: (item) => {
        if (!item.id) {
          return item.text
        }
        return item.text || item.id
      },
    })
  })
}

function buildURL() {
  var ext = `.${$("#fileTypeSelect").val()}`
  var slug = $("#proteinSlug").val()
  var title = $("#showName").prop("checked") ? "title=1" : ""
  var fill = $("#areaFill").prop("checked") ? "" : "fill=0"
  var transparent = $("#transCheck").prop("checked") ? "" : "transparent=0"
  var alpha = $("#opacitySlider").val() === "0.5" ? "" : `alpha=${$("#opacitySlider").val()}`
  var linewidth =
    $("#lineWidthSlider").val() === "1" ? "" : `linewidth=${$("#lineWidthSlider").val()}`
  var xlabels = $("#xAxis").prop("checked") ? "" : "xlabels=0"
  var ylabels = $("#yAxis").prop("checked") ? "ylabels=1" : ""
  var grid = $("#grid").prop("checked") ? "grid=1" : ""
  var xlim = ""
  if ($("#minXRange").val() !== "350" || $("#maxXRange").val() !== "750") {
    xlim = `xlim=${$("#minXRange").val()},${$("#maxXRange").val()}`
  }

  var newstring = `/spectra_img/${slug}${ext}`
  var argarray = [title, grid, xlabels, ylabels, fill, transparent, xlim, alpha, linewidth].filter(
    (d) => Boolean(d !== "")
  )
  if (argarray.length) {
    newstring += `?${argarray.join("&")}`
  }
  var fullurl = `${window.location.protocol}//${window.location.host}${newstring}`
  $("#activeImg").attr("src", newstring)
  $("#linktext").text(fullurl).attr("href", fullurl)
}

$("#spectra_url_form").submit((e) => {
  e.preventDefault()
})

$("#spectra_url_form input, #spectra_url_form select").change(function (_e) {
  if ($(this).hasClass("wave-range")) {
    if (!$(this).val()) {
      $(this).val($(this).hasClass("min-range") ? 350 : 750)
    }
    $(this).val(Math.min(Math.max($(this).val(), 200), 1600))
  }
  buildURL()
  if ($("#fileTypeSelect").val() === "pdf") {
    $("#activeImg").hide()
    $("#linktext").attr("download", true)
  } else {
    $("#activeImg").show()
    $("#linktext").removeAttr("download")
  }
})
if (document.getElementById("activeImg")) {
  buildURL()
}

// function get_color_group(exwave) {
//   if (exwave < 380) {
//     return "#C080FF"
//   } //
//   if (exwave < 420) {
//     return "#8080FF"
//   } // Blue
//   if (exwave < 473) {
//     return "#80FFFF"
//   } // Cyan
//   if (exwave < 507) {
//     return "#80FF80"
//   } // Green
//   if (exwave < 515) {
//     return "#CCFF80"
//   } // Yellow-Green
//   if (exwave < 531) {
//     return "#FFFF80"
//   } // Yellow
//   if (exwave < 555) {
//     return "#FFC080"
//   } // Orange
//   if (exwave < 600) {
//     return "#FFA080"
//   } // Red
//   if (exwave < 631) {
//     return "#FF8080"
//   } // Far red
//   if (exwave < 800) {
//     return "#B09090"
//   } // Near IR
// }

///////////////////////
// FORM BEHAVIOR
///////////////////////

// Hide unnecessary fiels in forms when the dark state checkbox is toggled
$(".dark_state_button input").click(function () {
  const neighbors = $(this).closest(".stateform_block").children(".hide_if_dark")
  if ($(this)[0].checked) {
    neighbors.hide()
    $(this).find("input[name*='max']").empty()
  } else {
    neighbors.show()
    $(this).closest(".stateform_block")
  }
})

$(() => {
  $(".dark_state_button input:checked").each(function () {
    const neighbors = $(this).closest(".stateform_block").children(".hide_if_dark")
    neighbors.hide()
  })
})

function reset_ipgid(hintstring) {
  if (hintstring) {
    $("#hint_id_ipg_id").html(hintstring)
  }
  $("#id_seq").prop("disabled", false)
  $("#hint_id_seq").html("Amino acid sequence (IPG ID is preferred)")
}

$("#id_ipg_id").change(function () {
  const ipg_id = $(this).val()
  const protein_uri =
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=json&rettype=fasta&id="
  const ipg_uri =
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=ipg&retmode=json&id="
  const fpbase_params = "&tool=fpbase&email=talley.lambert+fpbase@gmail.com"
  $.ajax({
    url: ipg_uri + ipg_id + fpbase_params,
    context: document.body,
    success: (data) => {
      if (!("result" in data)) {
        reset_ipgid(
          'NCBI <a href="https://www.ncbi.nlm.nih.gov/ipg/docs/about/">Identical Protein Group ID</a>'
        )
      } else if (ipg_id in data.result) {
        const accession = data.result[ipg_id].accession
        const title = data.result[ipg_id].title
        $("#hint_id_ipg_id").html(`IPG name: ${title}`)
        $.ajax({
          url: protein_uri + accession + fpbase_params,
          context: document.body,
          success: (data2) => {
            var lines = data2.split("\n")
            var seq = ""
            for (let i = 0; i < lines.length; i++) {
              if (lines[i].length !== 0 && lines[i][0] !== ">") {
                seq += lines[i]
              }
            }
            $("#id_seq").val(seq)
            //$("#id_seq").prop('disabled', true);
            //$("#hint_id_seq").html('Sequence input disabled when IPG ID provided')
          },
          error: (_data) => {
            reset_ipgid("Unrecognized IPG ID")
          },
        })
      }
    },
    error: (_data) => {
      reset_ipgid("Unrecognized IPG ID")
    },
  })
})

$("#proteinform #id_name").change(function () {
  var form = $(this).closest("form")
  $.ajax({
    method: "POST",
    url: form.data("validate-proteinname-url"),
    data: form.find("#id_slug:hidden, #id_name, [name='csrfmiddlewaretoken']").serialize(),
    dataType: "json",
    success: (data) => {
      if (data.is_taken) {
        const namelink = `<a href="${data.url}" style="text-decoration: underline;">${data.name}</a>`
        const message = `<strong>${namelink} already exists in the database.</strong>`
        $("#id_name").addClass("is-invalid")
        $("#div_id_name").addClass("has-danger")

        if ($("#error_1_id_name").length) {
          $("#error_1_id_name").html(message)
        } else {
          const span = $("<span/>", {
            id: "error_1_id_name",
            class: "invalid-feedback",
          }).append(message)
          $("#hint_id_name").before(span)
        }
      } else {
        if ($("#error_1_id_name").length) {
          $("#error_1_id_name").remove()
          $("#id_name").removeClass("is-invalid")
          $("#div_id_name").removeClass("has-danger")
        }
      }
    },
  })
})

$("#spectrum-submit-form #id_owner").change(function () {
  var form = $(this).closest("form")
  $.ajax({
    method: "POST",
    url: form.data("validate-owner-url"),
    data: form.find("#id_owner, [name='csrfmiddlewaretoken']").serialize(),
    dataType: "json",
    success: (data) => {
      if (data.similars.length) {
        let str = "<strong>Avoid duplicates.</strong> Similarly named existing spectra: "
        $.each(data.similars, (index, val) => {
          str = `${str}<span class="text-danger">${val.name}</span>`
          if (val.spectra.length) {
            str = `${str} (`
            $.each(val.spectra, (i, s) => {
              str = str + s
              if (i !== val.spectra.length - 1) {
                str = `${str}, `
              }
            })
            str = `${str})`
          }
          if (index !== data.similars.length - 1) {
            str = `${str}, `
          }
        })
        $("#hint_id_owner").html(str)
      } else {
        $("#hint_id_owner").html("Owner of the spectrum")
      }
    },
  })
})

// auto populate PMID after DOI input

$('input[id*="reference_doi"]').change(function () {
  const input = $(this)
  const small = input.parent().siblings('small[id*="reference_doi"]')
  const doi = input.val()
  var searchurl = "https://api.crossref.org/v1/works/"
  searchurl += doi

  // searchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="
  // searchurl += doi + "[doi]&retmode=json"

  $.ajax({
    url: searchurl,
    context: document.body,
    success: (data) => {
      if (data.status === "ok") {
        const year = data.message.issued["date-parts"]["0"]["0"]
        const author = data.message.author["0"].family
        const title = data.message.title[0].slice(0, 45)
        const citation = `${author} (${year}) ${title}...`
        small.html(citation)
      } else {
        small.html("DOI not found at Crossref")
      }
    },
    error: (_data) => {
      small.html("DOI not found at Crossref")
    },
  })
})

// // auto populate DOI after PMID input
// $('#id_reference_pmid').change(function(){
// 	pmid = $(this).val();
// 	console.log(pmid);
// 	searchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=";
// 	searchurl += pmid + "&retmode=json";
// 	$.ajax({
// 	    url: searchurl,
// 	    context: document.body,
// 	    success: function(data){
// 	    	console.log(data)
// 	    	if (pmid in data.result){
// 	    		articleids = data.result[pmid].articleids;
// 	    		doi = articleids.find(o => o.idtype === 'doi');
// 	    		if(doi){
// 	    			$('#id_reference_doi').val(doi.value);
// 	    		} else {
// 	    			$('#id_reference_doi').val('');
// 	    			$("#hint_id_reference_doi").html("DOI not found at Pubmed")
// 	    		}
// 	    		year = data.result[pmid].pubdate.split(" ")[0];
// 	    		author = data.result[pmid].sortfirstauthor.split(" ")[0];
// 	    		title = data.result[pmid].title.slice(0, 38);
// 	    		citation = author + " (" + year + ") " + title + "...";
// 	    		$("#hint_id_reference_doi").html(citation)
// 	    		$("#hint_id_reference_pmid").html("Pubmed ID located!")
// 	    	} else {
// 	    		$("#hint_id_reference_pmid").html("Unable to locate Pubmed ID")
// 	    		$('#id_reference_doi').val('');
// 	    		$("#hint_id_reference_doi").html("DOI not found at Pubmed")
// 	    	}
// 	    }
// 	});
// })

/////////////////// PROTEIN DETAIL PAGE //////////////////////

function chunkString(str, len) {
  const _size = Math.ceil(str.length / len)
  const _ret = new Array(_size)
  let _offset

  for (let _i = 0; _i < _size; _i++) {
    _offset = _i * len
    _ret[_i] = str.substring(_offset, _offset + len)
  }

  return _ret
}

function tooltipwrap(chunk, index, skipV2) {
  var out = ""
  for (let i = 0; i < chunk.length; i++) {
    let ind
    if (skipV2) {
      if (+index + i < 1) {
        ind = +index + i + 1
      } else if (+index + i === 1) {
        ind = "1a"
      } else if (+index + i > 1) {
        ind = +index + i
      }
    } else {
      ind = +index + i + 1
    }
    out +=
      '<span data-bs-toggle="tooltip" data-placement="top" title="' +
      chunk[i] +
      " " +
      ind +
      '">' +
      chunk[i] +
      "</span>"
  }
  return out
}

function formatAAseq(elem, breakpoint) {
  breakpoint = breakpoint || 10
  //clear any existing counts
  elem.find(".sequence_count").empty()
  // extract the string and chop it up into segments by breakpoint
  let skipV2
  if (elem.text().startsWith("MVSKGEEL")) {
    skipV2 = true
  }
  const words = chunkString(elem.text().replace(/ /g, ""), breakpoint)
  // clear the div
  elem.empty()
  //console.log(words)
  // create new divs for the formatted content
  const seqcount = $("<div class='sequence_count'></div>").appendTo(elem)
  const seqdiv = $("<div class='formatted_aminosquence'></div>").appendTo(elem)
  seqdiv.html(tooltipwrap(words[0], 0, skipV2))
  seqcount.append(`${1}<br>`)
  let height = seqdiv.height()
  for (let i = 1; i < words.length; i++) {
    const tippywords = tooltipwrap(words[i], i * 10, skipV2)
    seqdiv.html(`${seqdiv.html()} ${tippywords}`)
    if (seqdiv.height() > height) {
      // line break occured at iteration i
      //console.log(words[i])
      seqcount.append(`${i * breakpoint + 1}<br>`)
      height = seqdiv.height()
    }
  }
  elem.show()
  // Wait for Bootstrap to be available before initializing tooltips
  waitForBootstrap(() => {
    $('[data-bs-toggle="tooltip"]').tooltip({
      trigger: "hover",
      delay: { show: 200 },
    })
  })
}

$(() => {
  setTimeout(() => {
    // waiting is just a hack...
    $(".aminosequence").each(function () {
      formatAAseq($(this))
    })
  }, 1)
})

$(window).resize(() => {
  $(".aminosequence").each(function () {
    formatAAseq($(this))
  })
})

$("#refModalForm").submit(function (e) {
  var form = $(this).closest("form")
  $.ajax({
    type: "POST",
    url: form.attr("data-action-url"),
    data: form.serialize(),
    dataType: "json",
    success: (data) => {
      if (data.status === "success") {
        window.location.reload()
      }
    },
  })
  e.preventDefault() // avoid to execute the actual submit of the form.
})

$("#excerptModalForm")
  .off("submit")
  .on("submit", function (e) {
    e.preventDefault() // Prevent default FIRST
    var form = $(this).closest("form")
    var $submitBtn = form.find('button[type="submit"]')

    // Prevent double submission
    if ($submitBtn.prop("disabled")) return
    $submitBtn.prop("disabled", true)

    $.ajax({
      type: "POST",
      url: form.attr("data-action-url"),
      data: form.serialize(),
      dataType: "json",
      success: (data) => {
        if (data.status === "success") {
          window.location.reload()
        }
      },
      error: () => {
        $submitBtn.prop("disabled", false) // Re-enable on error
      },
    })
  })

function register_transition_form() {
  // Call formset plugin on the form elements inside the container
  $(".trans_formset_div .formset-form").formset({
    addText: "Add Transition",
    addCssClass: "btn btn-info mb-4",
    deleteCssClass: "transDelete",
    deleteText: icon("remove"),
    prefix: "transitions",
    formCssClass: "formset-form",
    processHidden: true,
  })
}

// This function is for showing the modal
$(() => {
  $("#show_transition_modal").click(function () {
    $.ajax({
      type: "GET",
      url: $(this).attr("data-action-url"),
      data: {},
      cache: false,
      success: (data, _status) => {
        $("#transitionForm").html(data)

        // Use Bootstrap's shown.bs.modal event to ensure modal is fully visible
        const $modal = $("#transitionModal")
        // Remove any previous handlers to avoid duplicates
        $modal.off("shown.bs.modal")
        // Register formset AFTER modal is fully shown
        $modal.one("shown.bs.modal", () => {
          register_transition_form()
        })
        // Now show the modal (which will trigger the event above)
        $modal.modal()
      },
      error: (_xhr, status, error) => {
        console.error("Transition form AJAX error:", status, error)
      },
    })
  })
})

$("#transitionForm").submit(function (e) {
  var form = $(this).closest("form")
  $.ajax({
    type: "POST",
    url: form.attr("data-action-url"),
    data: form.serialize(),
    cache: false,
    success: (_data, _status) => {
      window.location.reload()
    },
    error: (data, _status, _error) => {
      $("#transitionForm").html(data.responseText)
      register_transition_form()
    },
  })
  e.preventDefault() // avoid to execute the actual submit of the form.
})

$("#adminApprove, #adminRevert").submit(function (e) {
  var form = $(this).closest("form")
  $.ajax({
    type: "POST",
    url: form.attr("data-action-url"),
    data: form.serialize(),
    dataType: "json",
    success: (_data) => {
      window.location = form.data("success")
    },
  })
  e.preventDefault() // avoid to execute the actual submit of the form.
})

/////////////////// END PROTEIN DETAIL PAGE //////////////////////

////////////////// AJAX REMOVE FROM COLLECTION ////////////////////

$(document).ready(() => {
  $("a.object-flag").click(function (e) {
    e.preventDefault()
    var button = $(this)
    var wrapper = button.closest(".object-flag-wrapper")

    $.post({
      url: button.data("action-url"),
      data: {
        flagged: button.data("flagged"),
        target_model: button.data("model"),
        target_id: button.data("id"),
        csrfmiddlewaretoken: window.CSRF_TOKEN,
      },
      success: (response) => {
        if (response.status === "flagged") {
          button.data("flagged", 1)
          wrapper.addClass("is-flagged")
          button.data("original-title", "Admin: unflag")
        } else if (response.status === "unflagged") {
          button.data("flagged", 0)
          wrapper.removeClass("is-flagged")
          button.data("original-title", "Flag this excerpt for review")
        }
      },
    })
  })

  $(".btn.collection-remove").click(function (e) {
    var button = $(this)
    button.prop("disabled", true)
    $.ajax({
      url: button.attr("data-action-url"),
      type: "POST",
      data: {
        target_protein: button.data("object"),
        target_collection: button.data("collection"),
        csrfmiddlewaretoken: window.CSRF_TOKEN,
      },
      success: (response) => {
        if (response.status === "deleted") {
          button.closest("tr").remove()
        }
      },
    })
    e.preventDefault()
  })

  $(".collection-add-button").click(function (e) {
    //var button = $(this)
    $.ajax({
      type: "GET",
      url: $(this).attr("data-action-url"),
      dataType: "json",
      success: (data, _status) => {
        if ("members" in data) {
          const members = JSON.parse(data.members)
          if (members.length) {
            $("#currentmemberships").empty()
            $("<p>This protein is currently a member of these collections </p>").appendTo(
              "#currentmemberships"
            )
            const list = $("<ul>").appendTo("#currentmemberships")
            $.each(members, function (_e) {
              const li = $("<li>")
              $("<a>").html(this[0]).attr("href", this[1]).appendTo(li)
              li.appendTo(list)
            })
          }
        }
        if (!$("#collectionSelect").length) {
          // only retrieve once
          $("#collectionSelection").prepend(data.widget)
          $("#collectionModal").modal()
        }
      },
    })
    e.preventDefault()
  })

  $("#collectionForm").submit(function (e) {
    var form = $(this).closest("form")
    const data = form.serialize()
    $.ajax({
      type: "POST",
      url: form.attr("data-action-url"),
      data,
      cache: false,
      success: (_data, _status) => {},
    })
    $("#collectionModal").modal("hide")
    e.preventDefault() // avoid to execute the actual submit of the form.
  })
})

/////////////// ORGANISM BUTTON ///////////

$(() => {
  $("#id_parent_organism")
    .siblings(".input-group-append")
    .find(".select-add-button")
    .click(() => {
      $("#organismModal").modal()
    })

  $("#taxonomyModalForm").submit(function (e) {
    const form = $(this)
    const tax_id = form.find('input[name="taxonomy_id"]').val()
    const tax_uri =
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&retmode=json&id="
    $.ajax({
      url: tax_uri + tax_id,
      success: (data) => {
        if ("scientificname" in data.result[data.result.uids[0]]) {
          // successful fetch from NCBI
          const sci_name = data.result[data.result.uids[0]].scientificname
          $.ajax({
            type: "POST",
            url: form.attr("data-action-url"),
            data: form.serialize(),
            dataType: "json",
            success: (_data) => {
              $(`<option value="${tax_id}">${sci_name} </option>`).appendTo("#id_parent_organism")
              $(`#id_parent_organism option[value="${tax_id}"]`).prop("selected", true)
              $("#organismModal").modal("hide")
            },
          })
        } else {
          $("#id_taxonomy_id").addClass("is-invalid")
          $("#div_id_taxonomy_id").addClass("has-danger")
          $("#hint_id_taxonomy_id")
            .text("Could not find that Taxonomy id at NCBI!")
            .addClass("text-danger font-weight-bold")
            .removeClass("text-muted")
        }
      },
      error: (_data) => {},
    })
    e.preventDefault()
  })
})

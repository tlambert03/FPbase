// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/jquery-ajax-sentry.js" // Track jQuery AJAX errors

const $ = window.jQuery // jQuery loaded from CDN

// Helper to wait for Bootstrap plugins to be available
function _waitForBootstrap(callback) {
  if (typeof $.fn.tab !== "undefined") {
    callback()
  } else {
    // Retry after a short delay
    setTimeout(() => _waitForBootstrap(callback), 50)
  }
}

window.$ = window.jQuery = $
// select2 JS loaded from CDN in base.html - only import CSS
import "select2/dist/css/select2.css"
import "select2-bootstrap-5-theme/dist/select2-bootstrap-5-theme.min.css"
import "./js/jquery.formset"

// Mark this bundle for Sentry context
window.FPBASE = window.FPBASE || {}
window.FPBASE.currentBundle = "microscopeForm"

$("#microscopeform").submit((e) => {
  if (confirmError) {
    const multidichroics = $("select[id*='bs_filters']").filter(function () {
      return $(this).val().length > 1
    })
    const multiex = $("select[id*='ex_filters']").filter(function () {
      return $(this).val().length > 1
    })
    if (multidichroics.length > 0 || multiex.length > 0) {
      e.preventDefault()
      $("#sureModal")
        .one("hidden.bs.modal", () => {
          if ($(document.activeElement).attr("id") === "save") {
            confirmError = false
            $("#microscopeform").submit()
          }
        })
        .modal()
    }
  } else {
    $("#spinner").css("display", "inline")
  }
})

$(document).on("show.bs.modal", "#sureModal", () => {
  $("#sureModal").css("z-index", 1200)
})

$(document).on("hidden.bs.modal", "#sureModal", () => {
  $("#sureModal").css("z-index", 1039)
})

$("#id_detector, #id_light_source, #id_extra_cameras, #id_extra_lights").select2({
  theme: "bootstrap",
  containerCssClass: ":all:",
  placeholder: "---------",
  allowClear: true,
  width: "auto",
})

$(() => {
  $(".formset_div").formset({
    addText: "Add Optical Configuration",
    addCssClass: "btn btn-sm btn-info add-oc-button mb-1 mr-1",
    deleteCssClass: "oc-delete-button",
    deleteText: "&times;",
    processHidden: true,
    prefix: "optical_configs", // it's important that this match the string used for the
    // related_name attribute in the "Microscope" foreign key field
    // in the OpticalConfig model
  })

  if (formErrors) {
    // Wait for Bootstrap to be available before calling .tab()
    _waitForBootstrap(() => {
      $('#microscopeFormTabs a[href="#bulk"]').tab("show")
    })
  }
})

$("#chromaImportForm, #semrockImportForm").submit(function (e) {
  e.preventDefault() // avoid to execute the actual submit of the form.
  $("#footerSpinner").show()
  $("#footerFail").hide()
  $("#footerSuccess").hide()
  const form = $(this).closest("form")
  const brand = form.data("brand")
  $.ajax({
    type: "POST",
    url: form.attr("data-action-url"),
    data: form.serialize(),
    dataType: "json",
    success: (data) => {
      if (data.status) {
        const newdata = JSON.parse(data.spectra_options)
        $('.data-selector[data-category="f"]').append(
          $("<option>", {
            value: newdata.slug,
          }).text(newdata.name)
        )
        $(`#${brand}Input`).removeClass("is-invalid")
        $(`#${brand}Help`).removeClass("invalid-feedback").addClass("text-muted").text("Success!")
        $("#footerSpinner").hide()
        $("#footerFail").hide()
        $("#footerSuccess").show()
      } else {
        $(`#${brand}Input`).addClass("is-invalid")
        $(`#${brand}Help`)
          .removeClass("text-muted")
          .addClass("invalid-feedback")
          .text(`ERROR: ${data.message}`)
          .show()
        $("#footerSpinner").hide()
        $("#footerFail").show()
        $("#footerSuccess").hide()
      }
    },
  }).then((_d) => {
    $("#footerSpinner").hide()
    // $('#importModal').modal('hide')
  })
})

$(".importerClose").click(() => {
  $("#footerSpinner").hide()
  $("#footerFail").hide()
  $("#footerSuccess").hide()
})

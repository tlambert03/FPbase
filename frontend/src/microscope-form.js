// Vite modulepreload polyfill (must be first)
import "vite/modulepreload-polyfill"

// Initialize Sentry first to catch errors during module loading
import "./js/sentry-init.js"
import "./js/ajax-sentry.js" // Track jQuery AJAX errors
import { fetchWithSentry } from "./js/ajax-sentry"

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
  theme: "bootstrap-5",
  containerCssClass: ":all:",
  placeholder: "---------",
  allowClear: true,
  width: "resolve",
  selectionCssClass: "select2--small",
  dropdownCssClass: "select2--small",
})

$(() => {
  $(".formset_div").formset({
    addText: "Add Optical Configuration",
    addCssClass: "btn btn-sm btn-info add-oc-button mb-1 me-1",
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

  fetchWithSentry(form.attr("data-action-url"), {
    method: "POST",
    headers: {
      "Content-Type": "application/x-www-form-urlencoded",
      // Legacy header required by Django is_ajax() check in dual-purpose endpoints
      "X-Requested-With": "XMLHttpRequest",
    },
    body: form.serialize(),
  })
    .then((response) => response.json())
    .then((data) => {
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
    })
    .finally(() => {
      $("#footerSpinner").hide()
      // $('#importModal').modal('hide')
    })
})

$(".importerClose").click(() => {
  $("#footerSpinner").hide()
  $("#footerFail").hide()
  $("#footerSuccess").hide()
})

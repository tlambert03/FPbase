import $ from "jquery"
import "select2/dist/css/select2.css"
import "select2-theme-bootstrap4/dist/select2-bootstrap.css"
import "select2/dist/js/select2.full"
import "./js/jquery.formset"

$("#microscopeform").submit(function(e) {
  // eslint-disable-next-line
  if (confirmError) {
    const multidichroics = $("select[id*='bs_filters']").filter(function() {
      return $(this).val().length > 1
    })
    const multiex = $("select[id*='ex_filters']").filter(function() {
      return $(this).val().length > 1
    })
    if (multidichroics.length > 0 || multiex.length > 0) {
      e.preventDefault()
      $("#sureModal")
        .one("hidden.bs.modal", function() {
          if ($(document.activeElement).attr("id") === "save") {
            confirmError = false // eslint-disable-line
            $("#microscopeform").submit()
          }
        })
        .modal()
    }
  } else {
    $("#spinner").css("display", "inline")
  }
})

$(document).on("show.bs.modal", "#sureModal", function() {
  $("#sureModal").css("z-index", 1200)
})

$(document).on("hidden.bs.modal", "#sureModal", function() {
  $("#sureModal").css("z-index", 1039)
})

$(
  "#id_detector, #id_light_source, #id_extra_cameras, #id_extra_lights"
).select2({
  theme: "bootstrap",
  containerCssClass: ":all:",
  placeholder: "---------",
  allowClear: true,
  width: "auto",
})

$(function() {
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

  // eslint-disable-next-line
  if (formErrors) {
    $('#microscopeFormTabs a[href="#bulk"]').tab("show")
  }
})

$("#chromaImportForm, #semrockImportForm").submit(function(e) {
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
    success: function(data) {
      if (data.status) {
        const newdata = JSON.parse(data.spectra_options)
        $('.data-selector[data-category="f"]').append(
          $("<option>", {
            value: newdata.slug,
          }).text(newdata.name)
        )
        $(`#${brand}Input`).removeClass("is-invalid")
        $(`#${brand}Help`)
          .removeClass("invalid-feedback")
          .addClass("text-muted")
          .text("Success!")
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
  }).then(function(d) {
    $("#footerSpinner").hide()
    // $('#importModal').modal('hide')
  })
})

$(".importerClose").click(function() {
  $("#footerSpinner").hide()
  $("#footerFail").hide()
  $("#footerSuccess").hide()
})

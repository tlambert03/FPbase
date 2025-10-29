$(document).ready(() => {
  $("#add_remove_favorite").click(function (e) {
    var $obj = $(this)
    var _target_id = $obj.data("target").split("_")[1]
    $obj.prop("disabled", true)
    $.ajax({
      url: $obj.attr("data-action-url"),
      type: "POST",
      data: {
        target_model: $obj.data("model"),
        target_object_id: $obj.data("target").split("_")[1],
        csrfmiddlewaretoken: window.CSRF_TOKEN,
      },
      success: (response) => {
        if (response.status === "added") {
          $obj.find("i").removeClass("far").addClass("fas")
          $obj.find("span").text("Remove from favorites")
        } else if (response.status === "deleted") {
          $obj.find("i").removeClass("fas").addClass("far")
          $obj.find("span").text("Add to favorites")
        } else {
          console.warn("Unexpected response status:", response.status)
        }
        $obj.parent(".favit").children(".fav-count").text(response.fav_count)
        $obj.prop("disabled", false)
      },
      error: (_xhr, _status, error) => {
        console.error("Failed to toggle favorite:", error)
        $obj.prop("disabled", false)
      },
    })
    e.preventDefault() // avoid to execute the actual submit of the form.
  })

  $(".btn.unfave").click(function () {
    var $obj = $(this)
    $obj.prop("disabled", true)
    $.ajax({
      url: $obj.attr("data-action-url"),
      type: "POST",
      data: {
        target_model: $obj.data("model"),
        target_object_id: $obj.data("id"),
        csrfmiddlewaretoken: window.CSRF_TOKEN,
      },
      success: (response) => {
        if (response.status === "deleted") {
          $obj.parent().remove()
        }
      },
      complete: (_response) => {
        $obj.prop("disabled", false)
      },
    })
    e.preventDefault() // avoid to execute the actual submit of the form.
  })
})

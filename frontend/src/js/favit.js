$(document).ready(() => {
  $("#add_remove_favorite").click(function (e) {
    e.preventDefault()

    var $obj = $(this)
    var $favit = $obj.closest(".favit")
    $obj.prop("disabled", true)
    $.ajax({
      url: $obj.attr("data-action-url"),
      type: "POST",
      data: {
        target_model: $obj.data("model"),
        target_object_id: $obj.data("bsTarget").split("_")[1],
        csrfmiddlewaretoken: window.CSRF_TOKEN,
      },
      success: (response) => {
        if (response.status === "added") {
          $favit.addClass("is-favorited")
          $favit.attr("data-original-title", "Remove from favorites")
          $obj.find(".favorite-text").text("Remove from favorites")
        } else if (response.status === "deleted") {
          $favit.removeClass("is-favorited")
          $favit.attr("data-original-title", "Add to favorites")
          $obj.find(".favorite-text").text("Add to favorites")
        } else {
          console.warn("Unexpected response status:", response.status)
        }
        $favit.children(".fav-count").text(response.fav_count)
        $obj.prop("disabled", false)
      },
      error: (_xhr, _status, error) => {
        console.error("Failed to toggle favorite:", error)
        $obj.prop("disabled", false)
      },
    })
  })

  $(".btn.unfave").click(function (e) {
    e.preventDefault() // Prevent default FIRST
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
  })
})

$(document).ready(() => {
  $("#add_remove_favorite").click(function (e) {
    e.preventDefault() // avoid to execute the actual submit of the form.
    var $obj = $(this)
    var $favit = $obj.closest(".favit")
    $obj.prop("disabled", true)

    const formData = new URLSearchParams({
      target_model: $obj.data("model"),
      target_object_id: $obj.data("target").split("_")[1],
      csrfmiddlewaretoken: window.CSRF_TOKEN,
    })

    fetch($obj.attr("data-action-url"), {
      method: "POST",
      headers: {
        "Content-Type": "application/x-www-form-urlencoded",
      },
      body: formData,
    })
      .then((response) => response.json())
      .then((response) => {
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
      })
      .catch((error) => {
        console.error("Failed to toggle favorite:", error)
        $obj.prop("disabled", false)
      })
  })

  $(".btn.unfave").click(function (e) {
    e.preventDefault() // avoid to execute the actual submit of the form.
    var $obj = $(this)
    $obj.prop("disabled", true)

    const formData = new URLSearchParams({
      target_model: $obj.data("model"),
      target_object_id: $obj.data("id"),
      csrfmiddlewaretoken: window.CSRF_TOKEN,
    })

    fetch($obj.attr("data-action-url"), {
      method: "POST",
      headers: {
        "Content-Type": "application/x-www-form-urlencoded",
      },
      body: formData,
    })
      .then((response) => response.json())
      .then((response) => {
        if (response.status === "deleted") {
          $obj.parent().remove()
        }
      })
      .finally(() => {
        $obj.prop("disabled", false)
      })
  })
})

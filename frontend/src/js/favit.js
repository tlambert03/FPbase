import { farHeart, fasHeart } from "../icons/fa-icons.ts"

// Helper to create SVG icon element
function createHeartIcon(isFilled, className = "") {
  const icon = isFilled ? fasHeart : farHeart
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg")
  svg.setAttribute("viewBox", `0 0 ${icon.width} ${icon.height}`)
  svg.setAttribute("fill", "currentColor")
  svg.setAttribute("aria-hidden", "true")
  if (className) {
    svg.className.baseVal = className
  }

  const path = document.createElementNS("http://www.w3.org/2000/svg", "path")
  path.setAttribute("d", icon.path)
  svg.appendChild(path)

  return svg
}

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
          // Replace outline hearts with filled hearts
          $obj.find("svg.favorite-icon").each(function () {
            const classes = this.getAttribute("class")
            $(this).replaceWith(createHeartIcon(true, classes))
          })
          $obj.find("span").text("Remove from favorites")
        } else if (response.status === "deleted") {
          // Replace filled hearts with outline hearts
          $obj.find("svg.favorite-icon").each(function () {
            const classes = this.getAttribute("class")
            $(this).replaceWith(createHeartIcon(false, classes))
          })
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

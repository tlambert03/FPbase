// Wait for Bootstrap to be available (loaded from CDN)
function waitForBootstrap(callback) {
  if (typeof $.fn.tooltip !== "undefined" && typeof $.fn.popover !== "undefined") {
    callback()
  } else {
    // Retry after a short delay
    setTimeout(() => waitForBootstrap(callback), 50)
  }
}

$(() => {
  $("#protein-image .no-glow").fadeOut(2200)
  $("#protein-image .glow").fadeIn(2200)

  waitForBootstrap(() => {
    $('[data-toggle="tooltip"]').tooltip({
      trigger: "hover",
      delay: { show: 200 },
    })

    $('[data-toggle="popover"]').popover({ html: true })
    $(".popover-dismiss").popover({ trigger: "focus" })
  })
})

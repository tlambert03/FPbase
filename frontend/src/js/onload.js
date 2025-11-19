const _fade = (el, ms, out) => {
  if (out) {
    Object.assign(el.style, { transition: `opacity ${ms}ms`, opacity: "0" })
    setTimeout(() => {
      el.style.display = "none"
    }, ms)
  } else {
    el.style.display = ""
    el.style.opacity = "0"
    el.style.transition = `opacity ${ms}ms`
    el.offsetHeight // Force reflow
    el.style.opacity = "1"
  }
}

const waitForBootstrap = (cb) =>
  window.bootstrap ? cb() : setTimeout(() => waitForBootstrap(cb), 50)

const init = (Comp, sel, opts = {}) =>
  document.querySelectorAll(sel).forEach((el) => {
    new Comp(el, opts)
  })

document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("#protein-image .no-glow").forEach((el) => {
    _fade(el, 2200, true)
  })
  document.querySelectorAll("#protein-image .glow").forEach((el) => {
    _fade(el, 2200, false)
  })

  waitForBootstrap(() => {
    init(bootstrap.Tooltip, '[data-bs-toggle="tooltip"]', {
      trigger: "hover",
      delay: { show: 200 },
    })
    init(bootstrap.Popover, '[data-bs-toggle="popover"]', { html: true })
    init(bootstrap.Popover, ".popover-dismiss", { trigger: "focus" })
  })
})

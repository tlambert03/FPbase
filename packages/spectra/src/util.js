const debounce = (fn, time) => {
  let timeout

  return function (...args) {
    const functionCall = () => fn.apply(this, args)

    clearTimeout(timeout)
    timeout = setTimeout(functionCall, time)
  }
}

function isTouchDevice() {
  try {
    const prefixes = " -webkit- -moz- -o- -ms- ".split(" ")

    const mq = (query) => window.matchMedia(query).matches

    if (
      "ontouchstart" in window ||
      (typeof window.DocumentTouch !== "undefined" && document instanceof window.DocumentTouch)
    ) {
      return true
    }

    return mq(["(", prefixes.join("touch-enabled),("), "heartz", ")"].join(""))
  } catch (_e) {
    // console.error("(Touch detect failed)", e)
    return false
  }
}

export { debounce, isTouchDevice }

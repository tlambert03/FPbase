/**
 * Creates a debounced version of a function that delays execution until after
 * a specified time has elapsed since the last call.
 *
 * @param fn - The function to debounce
 * @param time - The delay in milliseconds
 * @returns A debounced version of the function
 */
export const debounce = <T extends (...args: never[]) => unknown>(
  fn: T,
  time: number
): ((...args: Parameters<T>) => void) => {
  let timeout: ReturnType<typeof setTimeout> | undefined

  return function (this: unknown, ...args: Parameters<T>) {
    const functionCall = () => fn.apply(this, args)

    if (timeout !== undefined) {
      clearTimeout(timeout)
    }
    timeout = setTimeout(functionCall, time)
  }
}

/**
 * Detects if the current device supports touch input.
 *
 * @returns true if the device supports touch, false otherwise
 */
export function isTouchDevice(): boolean {
  try {
    const prefixes = " -webkit- -moz- -o- -ms- ".split(" ")

    const mq = (query: string) => window.matchMedia(query).matches

    // Check for touch support using standard API
    if ("ontouchstart" in window) {
      return true
    }

    // Legacy check for old IE browsers (DocumentTouch is non-standard)
    const windowWithDocumentTouch = window as Window & {
      DocumentTouch?: { new (): unknown }
    }
    if (
      typeof windowWithDocumentTouch.DocumentTouch !== "undefined" &&
      windowWithDocumentTouch.DocumentTouch &&
      document instanceof windowWithDocumentTouch.DocumentTouch
    ) {
      return true
    }

    return mq(["(", prefixes.join("touch-enabled),("), "heartz", ")"].join(""))
  } catch (_e) {
    // console.error("(Touch detect failed)", e)
    return false
  }
}

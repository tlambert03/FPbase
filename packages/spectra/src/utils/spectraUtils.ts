import type { OwnerInfo } from "../types"

/**
 * Reshape raw spectra info from API into ownerInfo and spectraInfo
 * Used by App to populate metadata stores
 */
// biome-ignore lint/suspicious/noExplicitAny: Raw API data has dynamic structure
export function reshapeSpectraInfo(arr: any[]) {
  if (!arr) return { ownerInfo: {}, spectraInfo: {} }

  return arr.reduce(
    (prev, cur) => {
      if (!Object.hasOwn(prev.ownerInfo, cur.owner.slug)) {
        prev.ownerInfo[cur.owner.slug] = {
          category: cur.category.toUpperCase(),
          label: cur.owner.name,
          spectra: [],
          value: cur.owner.slug,
          url: cur.owner.url,
        }
      }
      prev.ownerInfo[cur.owner.slug].spectra.push({
        id: String(cur.id),
        subtype: cur.subtype.toUpperCase(),
      })
      prev.spectraInfo[cur.id] = {
        subtype: cur.subtype.toUpperCase(),
        owner: cur.owner.slug,
        label: cur.owner.name,
        category: cur.category.toUpperCase(),
      }
      return prev
    },
    // biome-ignore lint/suspicious/noExplicitAny: spectraInfo has dynamic structure from API
    { ownerInfo: {} as Record<string, OwnerInfo>, spectraInfo: {} as Record<string, any> }
  )
}

/**
 * Calculate area under curve using trapezoidal rule
 */
export function trapz(arr: [number, number][]): number {
  let sum = 0
  for (let i = 1; i < arr.length; i++) {
    sum += 0.5 * (arr[i][1] + arr[i - 1][1]) * (arr[i][0] - arr[i - 1][0])
  }
  return sum
}

/**
 * Calculate product of two spectra
 */
export function spectraProduct(
  ar1: [number, number][],
  ar2: [number, number][]
): [number, number][] {
  const output: [number, number][] = []
  const left = Math.max(ar1[0][0], ar2[0][0])
  const right = Math.min(ar1[ar1.length - 1][0], ar2[ar2.length - 1][0])

  const idx1 = ar1.findIndex((x) => x[0] >= left)
  const idx2 = ar2.findIndex((x) => x[0] >= left)

  for (let i = 0; i <= right - left; i++) {
    output.push([left + i, ar1[idx1 + i][1] * ar2[idx2 + i][1]])
  }
  return output
}

/**
 * Debounce a function
 */
// biome-ignore lint/suspicious/noExplicitAny: Generic function type requires any
export function debounce<T extends (...args: any[]) => any>(
  fn: T,
  time: number
): (...args: Parameters<T>) => void {
  let timeout: NodeJS.Timeout | null = null

  return (...args: Parameters<T>) => {
    if (timeout) clearTimeout(timeout)
    timeout = setTimeout(() => fn(...args), time)
  }
}

/**
 * Check if device is touch-enabled
 */
export function isTouchDevice(): boolean {
  try {
    const prefixes = " -webkit- -moz- -o- -ms- ".split(" ")
    const mq = (query: string) => window.matchMedia(query).matches

    if (
      "ontouchstart" in window ||
      (window.DocumentTouch && document instanceof window.DocumentTouch)
    ) {
      return true
    }

    return mq(["(", prefixes.join("touch-enabled),("), "heartz", ")"].join(""))
  } catch (_e) {
    return false
  }
}

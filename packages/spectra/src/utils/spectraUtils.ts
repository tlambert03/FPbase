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
    const curr = arr[i]!
    const prev = arr[i - 1]!
    sum += 0.5 * (curr[1] + prev[1]) * (curr[0] - prev[0])
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
  const left = Math.max(ar1[0]![0], ar2[0]![0])
  const right = Math.min(ar1[ar1.length - 1]![0], ar2[ar2.length - 1]![0])

  const idx1 = ar1.findIndex((x) => x[0] >= left)
  const idx2 = ar2.findIndex((x) => x[0] >= left)

  for (let i = 0; i <= right - left; i++) {
    output.push([left + i, ar1[idx1 + i]![1] * ar2[idx2 + i]![1]])
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

    if ("ontouchstart" in window) {
      return true
    }

    // Legacy DocumentTouch check (non-standard)
    const DocumentTouch = (window as typeof window & { DocumentTouch?: unknown }).DocumentTouch
    if (DocumentTouch && document instanceof (DocumentTouch as { new (): Document })) {
      return true
    }

    return mq(["(", prefixes.join("touch-enabled),("), "heartz", ")"].join(""))
  } catch (_e) {
    return false
  }
}

/**
 * Compute overlap spectrum from multiple input spectra
 * Shared logic used by both EfficiencyTable and useSpectraData
 */
export function computeOverlap(
  ...spectra: Array<{
    id: string
    customId?: string
    data: [number, number][]
    category: string
    owner: {
      id: string
      name: string
      qy?: number
      slug?: string
    }
  }>
): {
  id: string
  data: [number, number][]
  area: number
  category: string
  subtype: string
  color: string
  owner: {
    id: string
    name: string
    qy?: number
    slug?: string
  }
} {
  // Sort for consistent ID generation
  const sorted = [...spectra].sort((a, b) => {
    const aId = a.customId || a.id
    const bId = b.customId || b.id
    return String(aId).localeCompare(String(bId))
  })

  const idString = sorted.map((s) => s.customId || s.id).join("_")

  const ownerName = sorted.map((s) => s.owner.name).join(" & ")
  const ownerID = sorted.map((s) => s.owner.id).join("_")

  // Get qy from first spectrum that has it
  const qy = sorted.reduce<number | undefined>((acc, next) => next.owner.qy || acc, undefined)

  // Get slug from first protein/dye spectrum
  const slug = sorted.reduce<string | undefined>(
    (acc, next) => (["P", "D"].includes(next.category) ? next.owner.slug : acc),
    undefined
  )

  // Compute product of all spectra data
  const product = sorted.reduce<[number, number][]>(
    (acc, s) => (acc.length > 0 ? spectraProduct(acc, s.data) : s.data),
    []
  )

  return {
    id: idString,
    data: product,
    area: trapz(product),
    category: "O",
    subtype: "O",
    color: "#000000",
    owner: {
      id: ownerID,
      name: ownerName,
      qy,
      slug,
    },
  }
}

// In-memory autocomplete ranking for the site search box.
//
// score = text relevance + W * popularity
//
// Text relevance rewards (in decreasing order) an exact whole-field match
// ("egfp" -> EGFP), a whole-field prefix match, word-level exact/prefix
// matches, substring matches ("gfp" -> mEGFP), and typo matches. Popularity
// (0-1) is computed server-side from page views and favorites
// (see backend/proteins/search_index.py). Weights were tuned against real
// FPbase queries, so change them with care.

const DIACRITICS = /[\u0300-\u036f]/g
const SEPARATORS = /[^a-z0-9\u0370-\u03ff]+/

export const normalize = (s) =>
  String(s ?? "")
    .normalize("NFKD")
    .replace(DIACRITICS, "")
    .toLowerCase()

export const tokenize = (s) => normalize(s).split(SEPARATORS).filter(Boolean)

/** Optimal string alignment distance, returning `max + 1` as soon as it exceeds `max`. */
export function osa(a, b, max) {
  if (Math.abs(a.length - b.length) > max) return max + 1
  let prev2 = null
  let prev = Array.from({ length: b.length + 1 }, (_, j) => j)
  for (let i = 1; i <= a.length; i++) {
    const cur = [i]
    let rowMin = i
    for (let j = 1; j <= b.length; j++) {
      const cost = a[i - 1] === b[j - 1] ? 0 : 1
      let v = Math.min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost)
      if (i > 1 && j > 1 && a[i - 1] === b[j - 2] && a[i - 2] === b[j - 1]) {
        v = Math.min(v, prev2[j - 2] + 1)
      }
      cur.push(v)
      if (v < rowMin) rowMin = v
    }
    if (rowMin > max) return max + 1
    prev2 = prev
    prev = cur
  }
  return prev[b.length]
}

export const DEFAULTS = {
  quality: { exact: 1, prefix: 0.9, infix: 0.8, typo1: 0.45, typo2: 0.25 },
  tokenWeight: 2, // T: word-level match
  exactBonus: 2.5, // E: whole field equals the query
  prefixBonus: 0.3, // P: whole field starts with the query
  coverage: 1, // C: query length / field length
  popularity: 3.5, // W
  minInfix: 3,
  minTypo1: 3,
  minTypo2: 8,
  optionalWords: ["protein"],
  fallback: true, // if nothing matches every query word, allow partial matches
}

// match type of query token `t` against indexed word `w`
function matchWord(t, w, isId, o) {
  if (w === t) return "exact"
  if (w.startsWith(t)) return "prefix"
  if (isId) return null
  if (t.length >= o.minInfix && w.includes(t)) return "infix"
  const max = t.length >= o.minTypo2 ? 2 : t.length >= o.minTypo1 ? 1 : 0
  if (!max) return null
  let d = osa(t, w, max)
  // allow the typo to fall within a prefix of the word ("mchr" -> "mcherry")
  for (let L = t.length - 1; L <= t.length + 1 && d > 0; L++) {
    if (L > 0 && L < w.length) d = Math.min(d, osa(t, w.slice(0, L), max))
  }
  if (d > max) return null
  return d === 0 ? "prefix" : d === 1 ? "typo1" : "typo2"
}

export class SearchIndex {
  /**
   * @param {object[]} records - must have a numeric popularity `p` in [0, 1]
   * @param {{key: string, weight: number, id?: boolean}[]} fields - searchable
   *   fields; values may be strings or arrays of strings. `id` fields only
   *   match exactly or by prefix (no substrings or typos).
   * @param {object} [options] - overrides for `DEFAULTS`
   */
  constructor(records, fields, options = {}) {
    this.o = { ...DEFAULTS, ...options, quality: { ...DEFAULTS.quality, ...options.quality } }
    this.docs = records.map((record) => {
      const values = []
      for (const { key, weight, id = false } of fields) {
        const raw = record[key]
        for (const value of Array.isArray(raw) ? raw : [raw]) {
          if (value === undefined || value === null || value === "") continue
          const words = tokenize(value)
          values.push({ key, weight, id, value: String(value), words, compact: words.join("") })
        }
      }
      return { record, values, pop: record.p || 0 }
    })
  }

  /** @returns {{record: object, score: number}[]} best matches, highest score first */
  search(query, limit = 5) {
    let hits = this._search(query, false)
    if (!hits.length && this.o.fallback) hits = this._search(query, true)
    return hits.slice(0, limit)
  }

  _search(query, partial) {
    const o = this.o
    let tokens = tokenize(query)
    const qc = tokens.join("")
    if (!qc) return []
    const required = tokens.filter((t) => !o.optionalWords.includes(t))
    if (required.length) tokens = required
    const memos = tokens.map(() => new Map())
    const hits = []
    for (const doc of this.docs) {
      let tokenSum = 0
      let matched = 0
      tokens.forEach((t, ti) => {
        const memo = memos[ti]
        let best = 0
        for (const f of doc.values) {
          for (const w of f.words) {
            const key = f.id ? `#${w}` : w
            let m = memo.get(key)
            if (m === undefined) {
              m = matchWord(t, w, f.id, o)
              memo.set(key, m)
            }
            if (m) best = Math.max(best, o.quality[m] * f.weight)
          }
          // in multi-word queries, words may be glued together in the field ("neon green")
          if (!best && tokens.length > 1 && !f.id && f.compact.includes(t)) {
            best = o.quality.infix * f.weight
          }
        }
        if (best) matched++
        tokenSum += best
      })
      // whole-field matches use the query with separators removed ("td tomato")
      let fieldBonus = 0
      let coverage = 0
      for (const f of doc.values) {
        if (f.compact === qc) fieldBonus = Math.max(fieldBonus, o.exactBonus * f.weight)
        else if (f.compact.startsWith(qc))
          fieldBonus = Math.max(fieldBonus, o.prefixBonus * f.weight)
        if (!f.id && f.compact.includes(qc))
          coverage = Math.max(coverage, qc.length / f.compact.length)
      }
      const allMatched =
        matched === tokens.length || fieldBonus > 0 || (coverage > 0 && tokens.length > 1)
      // fallback for queries with no full match: at least half the words must match
      if (!allMatched && !(partial && matched >= tokens.length / 2)) continue
      let text = (o.tokenWeight * tokenSum) / tokens.length + fieldBonus + o.coverage * coverage
      if (!allMatched) text *= matched / tokens.length
      hits.push({ record: doc.record, score: text + o.popularity * doc.pop, doc })
    }
    hits.sort((a, b) => b.score - a.score)
    return hits
  }

  /** Values of non-`exclude` fields that match the query, e.g. to show "aka: sfGFP". */
  matchedValues(hit, query, exclude = []) {
    const tokens = tokenize(query)
    const qc = tokens.join("")
    const out = []
    for (const f of hit.doc.values) {
      if (exclude.includes(f.key)) continue
      const hitsWord = tokens.some((t) => f.words.some((w) => matchWord(t, w, f.id, this.o)))
      if (hitsWord || (qc && f.compact.includes(qc))) out.push({ key: f.key, value: f.value })
    }
    return out
  }
}

const HTML_ESCAPES = { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }
export const escapeHtml = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => HTML_ESCAPES[c])

/** HTML-escape `value`, wrapping the parts that match `query` in <em>. */
export function highlight(value, query, options = {}) {
  const o = { ...DEFAULTS, ...options }
  const text = String(value ?? "")
  // normalized text, with a map back to indices in the original string
  let norm = ""
  const map = []
  for (let i = 0; i < text.length; i++) {
    const n = normalize(text[i])
    for (let k = 0; k < n.length; k++) map.push(i)
    norm += n
  }
  map.push(text.length)
  const marked = new Array(text.length).fill(false)
  const mark = (start, end) => {
    for (let i = map[start]; i < map[end]; i++) marked[i] = true
  }
  const words = [...norm.matchAll(/[a-z0-9\u0370-\u03ff]+/g)]
  for (const t of tokenize(query)) {
    let found = false
    for (let i = norm.indexOf(t); i >= 0; i = norm.indexOf(t, i + 1)) {
      mark(i, i + t.length)
      found = true
    }
    if (found) continue
    // typo matches: mark the whole word
    for (const w of words) {
      if (matchWord(t, w[0], false, o)) mark(w.index, w.index + w[0].length)
    }
  }
  let html = ""
  for (let i = 0; i < text.length; i++) {
    if (marked[i] && !marked[i - 1]) html += "<em>"
    html += escapeHtml(text[i])
    if (marked[i] && !marked[i + 1]) html += "</em>"
  }
  return html
}

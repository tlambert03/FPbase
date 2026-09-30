// Run with: pnpm --filter fpbase test
import assert from "node:assert/strict"
import { test } from "node:test"
import { highlight, osa, SearchIndex, tokenize } from "./engine.js"

const FIELDS = [
  { key: "name", weight: 1 },
  { key: "aliases", weight: 0.85 },
  { key: "pdb", weight: 0.9, id: true },
]
const PROTEINS = [
  { name: "EGFP", aliases: ["enhanced GFP"], pdb: ["2Y0G"], p: 0.98 },
  { name: "avGFP", aliases: ["wtGFP", "GFP", "Green Fluorescent Protein"], p: 0.95 },
  { name: "mEGFP", p: 0.8 },
  { name: "EGFP203C", p: 0.1 },
  { name: "GFPxm163", p: 0 },
  { name: "mCherry", pdb: ["2H5Q"], p: 1 },
  { name: "mCherry2", p: 0.6 },
  { name: "mChartreuse", p: 0.3 },
  { name: "tdTomato", p: 0.9 },
  { name: "dTomato", p: 0.6 },
  { name: "mNeonGreen", p: 0.95 },
  { name: "Venus", p: 0.8 },
  { name: "mVenus", p: 0.85 },
  { name: "mScarlet-I", p: 0.8 },
]
const index = new SearchIndex(PROTEINS, FIELDS)
const names = (query, limit = 5) => index.search(query, limit).map((h) => h.record.name)

test("tokenize strips diacritics, case and separators", () => {
  assert.deepEqual(tokenize("mScarlet-I (Ñ2)"), ["mscarlet", "i", "n2"])
})

test("osa counts transpositions as one edit and bails early", () => {
  assert.equal(osa("egpf", "egfp", 2), 1)
  assert.equal(osa("abc", "abcdef", 1), 2)
})

test("prefixes of popular proteins rank first", () => {
  assert.equal(names("mch")[0], "mCherry")
  assert.equal(names("egf")[0], "EGFP")
})

test("popular substring matches beat obscure prefix matches", () => {
  const res = names("gfp")
  assert.deepEqual(res.slice(0, 2), ["avGFP", "EGFP"])
  assert.ok(res.indexOf("mEGFP") < res.indexOf("GFPxm163"))
  assert.ok(names("egfp").indexOf("mEGFP") < names("egfp").indexOf("EGFP203C"))
})

test("exact whole-name matches beat more popular partial matches", () => {
  assert.equal(names("venus")[0], "Venus")
  assert.equal(names("mcherry")[0], "mCherry")
})

test("typos are tolerated", () => {
  assert.equal(names("mchery")[0], "mCherry")
  assert.equal(names("egpf")[0], "EGFP")
  assert.equal(names("tdtomatto")[0], "tdTomato")
})

test("separators in the query or the name don't matter", () => {
  assert.equal(names("td tomato")[0], "tdTomato")
  assert.equal(names("neon green")[0], "mNeonGreen")
  assert.equal(names("mscarleti")[0], "mScarlet-I")
  assert.equal(names("mscarlet i")[0], "mScarlet-I")
})

test("the word 'protein' is optional", () => {
  assert.equal(names("green fluorescent protein")[0], "avGFP")
})

test("ids match exactly or by prefix only", () => {
  assert.equal(names("2h5q")[0], "mCherry")
  assert.equal(names("2y0")[0], "EGFP")
  assert.deepEqual(names("h5q"), [])
})

test("no match returns nothing, partial matches need half the words", () => {
  assert.deepEqual(names("zzzzzz"), [])
  assert.deepEqual(names("alexa fluor 488"), [])
  assert.deepEqual(names("monomeric venus").sort(), ["Venus", "mVenus"])
})

test("matchedValues reports the non-name fields that matched", () => {
  const [hit] = index.search("wtgfp", 1)
  assert.deepEqual(index.matchedValues(hit, ["name"]), [{ key: "aliases", value: "wtGFP" }])
})

test("highlight escapes html and marks matches", () => {
  assert.equal(highlight("mCherry", "cher"), "m<em>Cher</em>ry")
  assert.equal(highlight("td Tomato", "td tomato"), "<em>td</em> <em>Tomato</em>")
  assert.equal(highlight("<b>x</b>", "b"), "&lt;<em>b</em>&gt;x&lt;/<em>b</em>&gt;")
  assert.equal(highlight("Ñeon", "neon"), "<em>Ñeon</em>")
  assert.equal(highlight("mCherry", "mchery"), "<em>mCherry</em>")
})

const DYES = new SearchIndex(
  [
    { name: "Alexa Fluor 488", p: 1 },
    { name: "Alexa Fluor 647", p: 0.95 },
    { name: "Alexa Fluor 430", p: 0.3 },
    { name: "CF647", p: 0.4 },
    { name: "Brilliant Violet 421", p: 0.5 },
  ],
  [{ key: "name", weight: 1 }],
  { synonyms: { af: "alexa fluor", bv: "brilliant violet" } }
)
const dyes = (query) => DYES.search(query, 3).map((h) => h.record.name)

test("letters glued to digits are split when nothing else matches", () => {
  assert.equal(dyes("alexa488")[0], "Alexa Fluor 488")
  assert.deepEqual(names("qqq123"), [])
})

test("synonyms expand abbreviations, alone or before digits", () => {
  assert.equal(dyes("af647")[0], "Alexa Fluor 647")
  assert.equal(dyes("af 488")[0], "Alexa Fluor 488")
  assert.equal(dyes("bv421")[0], "Brilliant Violet 421")
  assert.equal(dyes("cf647")[0], "CF647")
})

test("hits carry the expanded tokens for highlighting", () => {
  const [hit] = DYES.search("af647", 1)
  assert.equal(highlight(hit.record.name, hit.tokens), "<em>Alexa</em> <em>Fluor</em> <em>647</em>")
})

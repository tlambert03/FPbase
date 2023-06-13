import React from "react"

const zip = (arr, ...arrs) => {
  return arr.map((val, i) => arrs.reduce((a, arr) => [...a, arr[i]], [val]))
}

const chunkString = (str, length) =>
  str.match(new RegExp(".{1," + length + "}", "g"))

function BlastStatsTable({ hit }) {
  return (
    <table
      className={"small border-top border-bottom mb-2 flip-scroll"}
      style={{ width: "100%" }}
    >
      <tbody>
        <tr className="text-muted">
          <th>Score</th>
          <th>Expect</th>
          <th>Identities</th>
          <th>Positives</th>
          <th>Gaps</th>
          {hit.query_frame ? <th>Frame</th> : null}
        </tr>
        <tr>
          <td>
            {hit.bit_score} bits({hit.score})
          </td>
          <td>{hit.evalue.toExponential(2)}</td>
          <td>
            {hit.identity}/{hit.align_len}(
            {Math.round((1000 * hit.identity) / hit.align_len) / 10}
            %)
          </td>
          <td>
            {hit.positive}/{hit.align_len}(
            {Math.round((1000 * hit.positive) / hit.align_len) / 10}
            %)
          </td>
          <td>
            {hit.gaps}/{hit.align_len}(
            {Math.round((1000 * hit.gaps) / hit.align_len) / 10}%)
          </td>
          {hit.query_frame ? <td>{hit.query_frame}</td> : null}
        </tr>
      </tbody>
    </table>
  )
}

function FormattedBlastAlignment({ hit, lineWidth }) {
  const text = zip(
    chunkString(hit.qseq, lineWidth),
    chunkString(hit.midline, lineWidth),
    chunkString(hit.hseq, lineWidth)
  ).map(([q, m, h], index) => {
    if (q !== h) {
      let qsplit = q.split("")
      let hsplit = h.split("")
      let msplit = m.split("")
      qsplit.forEach((letter, index) => {
        if (letter.toUpperCase() !== hsplit[index].toUpperCase()) {
          qsplit[index] = `<span class="mismatch">${qsplit[index]}</span>`
          hsplit[index] = `<span class="mismatch">${hsplit[index]}</span>`
          msplit[index] = `<span class="mismatch">${msplit[index]}</span>`
        }
      })
      q = qsplit.join("")
      h = hsplit.join("")
      m = msplit.join("")
    }
    return [
      `Query   ${String(index * lineWidth + hit.query_from).padEnd(5)}${q}`,
      `             ${m}`,
      `Sbjct   ${String(index * lineWidth + hit.hit_from).padEnd(5)}${h}`,
      "" // adds a space between triplets
    ].join("\n")
  })
  return <pre dangerouslySetInnerHTML={{ __html: text.join("\n") }} />
}

function BlastHitSummary({ data }) {
  const accession = data.description[0].accession
  const title = data.description[0].title
  const id = data.description[0].id

  return (
    <div id={"dln_" + accession} className="dlfRow mt-4">
      <h5 style={{ fontWeight: "bold", color: "#5b616b" }}>{title}</h5>
      <div className="small">
        <label className="mr-1 font-weight-bold text-muted">FPbase ID:</label>
        <a
          href={`/protein/${accession}`}
          target="_blank"
          title={`Go to ${title} at FPbase`}
        >
          {accession}
        </a>
        <span className="ml-4">
          <label className="mr-1 font-weight-bold text-muted">Length: </label>
          {data.len}
        </span>
      </div>
    </div>
  )
}

function BlastHit({ data }) {
  return (
    <div className={"mt-2 border-bottom"}>
      <BlastHitSummary data={data} />
      <BlastStatsTable hit={data.hsps[0]} />
      <FormattedBlastAlignment hit={data.hsps[0]} lineWidth={60} />
    </div>
  )
}

function BlastReportAlignments({ report }) {
  const hits = report.search.hits

  return (
    <div className="mx-2">
      {hits.map(hit => (
        <BlastHit key={hit.num} data={hit} />
      ))}
    </div>
  )
}

export default BlastReportAlignments

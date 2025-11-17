import { Box, FormControl, Grid, InputLabel, MenuItem, Select, Typography } from "@mui/material"
import { useState } from "react"
import { fetchWithSentry } from "../../../frontend/src/js/ajax-sentry"

const $ = window.jQuery

import InputForm from "./InputForm.jsx"
import BlastReport from "./ReportView.jsx"

function ReportSelect({ reports, binary, index, onChange }) {
  const unit = binary === "blastp" ? "aa" : "nt"

  function handleChange(event) {
    onChange(event.target.value)
  }

  return (
    <Box sx={{ mt: 4 }}>
      <Grid container spacing={2} alignItems="center">
        <Grid item xs={12} sm={3} md={2}>
          <Typography variant="body1" sx={{ fontWeight: "bold", color: "#5b616b" }}>
            Results for:
          </Typography>
        </Grid>
        <Grid item xs={12} sm={9} md={10}>
          <FormControl fullWidth>
            <InputLabel id="report-select-label">Select Report</InputLabel>
            <Select
              labelId="report-select-label"
              id="report-select"
              value={index}
              label="Select Report"
              onChange={handleChange}
            >
              {reports.map((item, i) => (
                <MenuItem key={item.report.results.search.query_id} value={i}>
                  {`${item.report.results.search.query_id}: ${
                    item.report.results.search.query_title
                  } (${item.report.results.search.query_len}${unit})`}
                </MenuItem>
              ))}
            </Select>
          </FormControl>
        </Grid>
      </Grid>
    </Box>
  )
}

function App() {
  const [results, setResults] = useState([])
  const [binary, setBinary] = useState("blastp")
  const [reportIndex, setReportIndex] = useState(0)
  const notDNA = new RegExp(/[^(A|T|C|G)]/g)

  function handleSubmit(e) {
    e.preventDefault()
    setReportIndex(0)

    const seqLetters = $(e.target)[0][1]
      .value.toUpperCase()
      .replace(/^>.*$/gm, "")
      .replace(/(?:\r\n|\r|\n)/g, "")

    const bin = notDNA.test(seqLetters) ? "blastp" : "blastx"
    setBinary(bin)

    fetchWithSentry("", {
      method: "POST",
      headers: {
        "Content-Type": "application/x-www-form-urlencoded",
        // Legacy header required by Django is_ajax() check in dual-purpose endpoints
        "X-Requested-With": "XMLHttpRequest",
      },
      body: `${$(e.target).serialize()}&binary=${bin}`,
    })
      .then((response) => response.json())
      .then((data) => {
        if (data.status === 200) {
          setResults(data.blastResult)
        } else if (data.status === 500) {
          console.error(data.error)
          alert(
            "There was an error processing your input.  Please double check that it is an amino acid or nucleotide sequence, or multiple sequences in FASTA format"
          )
        }
      })
  }

  function handleChangeReport(index) {
    setReportIndex(index)
  }

  const initialValue = new URLSearchParams(window.location.search).get("query")

  return (
    <div>
      <InputForm handleSubmit={handleSubmit} initialValue={initialValue} />
      {results.length > 1 ? (
        <ReportSelect
          reports={results}
          binary={binary}
          index={reportIndex}
          onChange={handleChangeReport}
        />
      ) : null}
      {results.length ? <BlastReport report={results[reportIndex]} /> : null}
      {results.length ? (
        <div className="mt-4 text-muted">
          <p className="text-center small">
            Version: {results[0].report.version}
            <br />
            FPbase Sequence Database ({results[0].report.results.search.stat.db_num} sequences,{" "}
            {results[0].report.results.search.stat.db_len.toLocaleString("en")} total letters)
            <br />
            Matrix: {results[0].report.params.matrix}
          </p>
          <p className="small text-muted mt-4">
            <span className="font-weight-bold">Reference:</span>
            <br />
            Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng
            Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new
            generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.
          </p>
          <p className="small">
            <span className="font-weight-bold">
              Reference for compositional score matrix adjustment:
            </span>
            <br />
            Stephen F. Altschul, John C. Wootton, E. Michael Gertz, Richa Agarwala, Aleksandr
            Morgulis, Alejandro A. Schäffer, and Yi-Kuo Yu (2005) "Protein database searches using
            compositionally adjusted substitution matrices", FEBS J. 272:5101-5109.
          </p>
        </div>
      ) : null}
    </div>
  )
}

export default App

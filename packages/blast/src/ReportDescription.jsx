import React from "react"
import Table from "@material-ui/core/Table"
import TableBody from "@material-ui/core/TableBody"
import TableCell from "@material-ui/core/TableCell"
import TableHead from "@material-ui/core/TableHead"
import TableRow from "@material-ui/core/TableRow"
import Paper from "@material-ui/core/Paper"
import { makeStyles } from "@material-ui/core/styles"

function fpbaseLink(accession) {
  return (
    <a href={`/protein/${accession}`} target="_blank" rel="noopener noreferrer">
      {accession}
    </a>
  )
}

const useStyles = makeStyles(theme => ({
  table: {
    marginTop: "10px",
    minWidth: 650
  }
}))

function BlastReportDescription({ report, onClick }) {
  const rows = report.search.hits.map((elem, i) => {
    return {
      key: elem.num,
      title: elem.description[0].title,
      accession: elem.description[0].accession,
      evalue: elem.hsps[0].evalue.toExponential(2),
      bit_score: elem.hsps[0].bit_score,
      len: elem.len,
      gaps: elem.hsps[0].gaps,
      identity: elem.hsps[0].identity,
      align_len: elem.hsps[0].align_len
    }
  })

  const classes = useStyles()
  return (
    <Paper style={{ overflowX: "scroll" }}>
      <Table className={classes.table} size="small">
        <TableHead>
          <TableRow>
            <TableCell align="left">Description</TableCell>
            <TableCell align="right">Score&nbsp;(Bits)</TableCell>
            <TableCell align="right">E&nbsp;Value</TableCell>
            <TableCell align="right">Identities&nbsp;(%)</TableCell>
            <TableCell align="right">Gaps</TableCell>
            <TableCell align="right">Accession</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {rows.map(row => (
            <TableRow key={row.key}>
              <TableCell align="left">
                <a
                  className="text-secondary"
                  href={row.accession}
                  onClick={onClick}
                  title={`show alignment of query with ${row.title}`}
                >
                  {row.title}
                </a>
              </TableCell>
              <TableCell align="right">{row.bit_score}</TableCell>
              <TableCell align="right">{row.evalue}</TableCell>
              <TableCell align="right">
                {row.identity}/{row.align_len} (
                {Math.round((1000 * row.identity) / row.align_len) / 10}
                %)
              </TableCell>
              <TableCell align="right">{row.gaps}</TableCell>
              <TableCell align="right">{fpbaseLink(row.accession)}</TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </Paper>
  )
}

export default BlastReportDescription

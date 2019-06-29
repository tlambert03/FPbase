import React from "react"
import Table from "@material-ui/core/Table"
import TableBody from "@material-ui/core/TableBody"
import TableCell from "@material-ui/core/TableCell"
import TableHead from "@material-ui/core/TableHead"
import TableRow from "@material-ui/core/TableRow"
import Paper from "@material-ui/core/Paper"
import { makeStyles } from "@material-ui/core/styles"
import useSpectralData from "./Components/useSpectraData"


const useStyles = makeStyles(theme => ({
  table: {
    marginTop: "10px",
    minWidth: 650
  }
}))

function trapz(arr) {
  // approximate area under curve as series of trapezoids
  // arr = [[wave, data], ...]
  let sum = 0
  for (let i = 1; i < arr.length; i++) {
    sum += 0.5 * (arr[i][1] + arr[i - 1][1]) * (arr[i][0] - arr[i - 1][0])
  }
  return sum
}

function spectraProduct(ar1, ar2) {
  // calculate product of two spectra.values
  // these assume monotonic increase w/ step = 1
  const output = []
  const left = Math.max(ar1[0][0], ar2[0][0]) // the min wavelength shared by both arrays
  const right = Math.min(ar1[ar1.length - 1][0], ar2[ar2.length - 1][0]) // the max wavelength shared by both arrays
  if (left >= right) return []

  const offsetA1 = left - ar1[0][0]
  const offsetA2 = left - ar2[0][0]
  for (let i = 0; i < right - left; i++) {
    output.push([left + i, ar1[offsetA1 + i][1] * ar2[offsetA2 + i][1]])
  }
  return output
}

const EfficiencyTable = ({ activeSpectra, spectraInfo }) => {
  const data = useSpectralData(activeSpectra)

  const classes = useStyles()
  const filters = data.filter(({ category }) => category === "F")
  const emSpectra = data.filter(({ subtype }) => subtype === "EM")

  const rows = []
  const headers = []
  emSpectra.forEach(({ owner }) => {
    headers.push(owner.name)
  })
  filters.forEach(filter => {
    const row = [{ key: filter.id, data: filter.owner.name }]
    emSpectra.forEach(emSpectra => {
      row.push({
        data: (
          trapz(spectraProduct(filter.data, emSpectra.data)) / emSpectra.area
        ).toFixed(2),
        key: `${filter.id}_${emSpectra.id}`
      })
    })
    rows.push(row)
  })
  console.log(rows)
  return (
    <Paper style={{ overflowX: "scroll", marginBottom: 200 }}>
      {/* <div>filters: {JSON.stringify(data)}</div> */}
      <Table className={classes.table} size="small">
        <TableHead>
          <TableRow>
            <TableCell />
            {headers.map(header => (
              <TableCell key={header}>{header}</TableCell>
            ))}
          </TableRow>
        </TableHead>
        <TableBody>
          {rows.map(row => (
            <TableRow key={row[0].key}>
              {row.map(({ data, key }, index) =>
                index === 0 ? (
                  <TableCell key={key} align="right">
                    {data}
                  </TableCell>
                ) : (
                  <TableCell key={key} align="left">
                    {data}
                  </TableCell>
                )
              )}
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </Paper>
  )
}

export default EfficiencyTable
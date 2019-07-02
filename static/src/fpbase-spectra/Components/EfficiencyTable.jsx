import React, { useState } from "react"
import { makeStyles } from "@material-ui/core/styles"
import useSpectralData from "./useSpectraData"
import { trapz } from "../util"
import { Typography } from "@material-ui/core"
import MaterialTable from "material-table"
import Export from "@material-ui/icons/SaveAlt"
import FirstPage from "@material-ui/icons/FirstPage"
import PreviousPage from "@material-ui/icons/ChevronLeft"
import NextPage from "@material-ui/icons/ChevronRight"
import LastPage from "@material-ui/icons/LastPage"
import Search from "@material-ui/icons/Search"
import ResetSearch from "@material-ui/icons/Clear"
import Filter from "@material-ui/icons/FilterList"
import Shuffle from "@material-ui/icons/Shuffle"

const TABLE_ICONS = {
  Export,
  FirstPage,
  PreviousPage,
  NextPage,
  LastPage,
  Search,
  ResetSearch,
  Filter
}

const TABLE_OPTIONS = {
  search: false,
  paging: false,
  exportButton: true,
  exportAllData: true,
  padding: "dense"
}

const useStyles = makeStyles(theme => ({
  table: {
    marginTop: "10px",
    minWidth: 650
  },
  description: {
    padding: "8px 20px"
  }
}))

const OverlapCache = {}

function getOverlap(...args) {
  const idString = args
    .map(arg => arg.customId || arg.id)
    .sort((a, b) => a - b)
    .join("_")
  if (!(idString in OverlapCache)) {
    const product = spectraProduct(...args.map(({ data }) => data))
    OverlapCache[idString] = {
      data: product,
      area: trapz(product),
      key: idString
    }
  }
  return OverlapCache[idString]
}

function spectraProduct(ar1, ar2) {
  console.log("PRODUCT")
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

const EfficiencyTable = ({ activeSpectra, initialTranspose }) => {
  const [transposed, setTransposed] = useState(initialTranspose)
  const data = useSpectralData(activeSpectra)

  const classes = useStyles()
  const filters = data.filter(
    ({ category, subtype }) => category === "F" && subtype !== "BX"
  )
  const emSpectra = data.filter(({ subtype }) => subtype === "EM")
  
  console.log(filters)
  let rows, headers

  const cellStyle = { fontSize: "1rem" }
  if (transposed) {
    rows = []
    headers = [{ title: "Fluorophore", field: "fluor" }]
    filters.forEach(({ owner }) => {
      headers.push({
        title: owner.name,
        field: owner.id,
        type: "numeric",
        cellStyle
      })
    })
    emSpectra.forEach(fluor => {
      const row = { fluor: fluor.owner.name }
      filters.forEach(filter => {
        const overlap = getOverlap(fluor, filter)
        row[filter.owner.id] = (100 * overlap.area / fluor.area).toFixed(1)
      })
      rows.push(row)
    })
  } else {
    rows = []
    headers = [{ title: "Filter", field: "filter" }]
    emSpectra.forEach(({ owner }) => {
      headers.push({
        title: owner.name,
        field: owner.id,
        type: "numeric",
        cellStyle
      })
    })
    filters.forEach(filter => {
      const row = { filter: filter.owner.name }
      emSpectra.forEach(emSpectra => {
        const overlap = getOverlap(filter, emSpectra)
        row[emSpectra.owner.id] = (100 * overlap.area / emSpectra.area).toFixed(1)
      })
      rows.push(row)
    })
  }

  if (rows.length < 1 || headers.length < 2) {
    return (
      <div className={classes.description}>
        <Typography variant="h6">Efficiency Table</Typography>
        <Typography variant="body1">
          Add at least on filter and one fluorophore, and this tab will show a
          table of collection efficiency (sometimes called "spillover") for each
          filter/fluorophore combination
        </Typography>
      </div>
    )
  }

  return (
    <div className="efficiency-table">
      <MaterialTable
        columns={headers}
        data={rows}
        options={TABLE_OPTIONS}
        icons={TABLE_ICONS}
        title="Collection Efficiency (%)"
        actions={[
          {
            icon: Shuffle,
            tooltip: "Transpose",
            isFreeAction: true,
            onClick: event => setTransposed(transposed => !transposed)
          }
        ]}
      />
    </div>
  )
}

export default EfficiencyTable

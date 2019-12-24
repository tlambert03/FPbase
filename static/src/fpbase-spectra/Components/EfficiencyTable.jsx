import React, { useState, useEffect } from "react"
import { makeStyles } from "@material-ui/core/styles"
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
import { useMutation, useQuery, useApolloClient } from "@apollo/react-hooks"
import Button from "@material-ui/core/Button"
import gql from "graphql-tag"
import { trapz } from "../util"
import useSpectralData from "./useSpectraData"
import { GET_ACTIVE_OVERLAPS } from "../client/queries"

class ErrorBoundary extends React.Component {
  static getDerivedStateFromError(error) {
    return { hasError: true }
  }

  componentDidCatch(error, info) {}

  render() {
    const { children } = this.props
    return children
  }
}

const TABLE_ICONS = {
  Export,
  FirstPage,
  PreviousPage,
  NextPage,
  LastPage,
  Search,
  ResetSearch,
  Filter,
}

const TABLE_OPTIONS = {
  search: false,
  paging: false,
  exportButton: true,
  exportAllData: true,
  padding: "dense",
}

const useStyles = makeStyles(theme => ({
  table: {
    marginTop: "10px",
    minWidth: 650,
  },
  description: {
    padding: "8px 20px",
  },
}))

window.OverlapCache = {}

function numStringSort(a, b) {
  if (Number.isNaN(parseFloat(a))) {
    if (Number.isNaN(parseFloat(b))) {
      return a - b
    }
    return -1
  }
  return a - b
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

function getOverlap(...args) {
  const idString = args
    .map(arg => arg.customId || arg.id)
    .sort(numStringSort)
    .join("_")

  const ownerName = args.map(({ owner }) => owner.name).join(" & ")
  const ownerID = args
    .map(({ owner }) => owner.id)
    .sort(numStringSort)
    .join("_")
  const qy = args.reduce((acc, next) => next.owner.qy || acc, null)
  const slug = args.reduce(
    (acc, next) => (["P", "D"].includes(next.category) ? next.owner.slug : acc),
    null
  )

  if (!(idString in window.OverlapCache)) {
    const product = spectraProduct(...args.map(({ data }) => data))
    window.OverlapCache[idString] = {
      data: product,
      area: trapz(product),
      id: idString,
      category: "O",
      subtype: "O",
      color: "#000000",
      owner: { id: ownerID, name: ownerName, qy, slug },
    }
  }
  return window.OverlapCache[idString]
}

const EfficiencyTable = ({ initialTranspose }) => {
  const [transposed, setTransposed] = useState(initialTranspose)
  const [[rows, headers], setTableData] = useState([[], []])
  const classes = useStyles()
  const spectraData = useSpectralData()
  const client = useApolloClient()

  useEffect(() => {
    return () => {
      client.writeData({ data: { activeOverlaps: [] } })
    }
  }, [client])

  useEffect(() => {
    async function updateTableData() {
      const filters = spectraData.filter(
        ({ category, subtype }) => category === "F" && subtype !== "BX"
      )
      const emSpectra = spectraData.filter(({ subtype }) => subtype === "EM")

      // untransposed columns represent different fluors
      let colItems = emSpectra
      let rowItems = filters
      let newHeaders = [{ title: "Filter", field: "field" }]
      const newRows = []
      if (transposed) {
        // columns represent different filters
        colItems = filters
        rowItems = emSpectra
        newHeaders = [{ title: "Fluorophore", field: "field" }]
      }

      colItems.forEach(({ owner }) => {
        newHeaders.push({
          title: owner.name,
          field: owner.id,
          type: "numeric",
          cellStyle: { fontSize: "1rem" },
          render: rowData => {
            const overlapID = rowData[`${owner.id}_overlapID`]
            return (
              <OverlapToggle id={overlapID}>{rowData[owner.id]}</OverlapToggle>
            )
          },
        })
      })

      rowItems.forEach(rowItem => {
        const row = { field: rowItem.owner.name }
        colItems.forEach(colItem => {
          const overlap = getOverlap(rowItem, colItem)
          const fluor = transposed ? rowItem : colItem
          row[colItem.owner.id] = ((100 * overlap.area) / fluor.area).toFixed(1)
          row[`${colItem.owner.id}_overlapID`] = overlap.id
        })
        newRows.push(row)
      })
      setTableData([newRows, newHeaders])
    }

    updateTableData()
  }, [client, spectraData, transposed])

  if (rows.length < 1 || headers.length < 2) {
    return (
      <div className={classes.description}>
        <Typography variant="h6">Efficiency Table</Typography>
        <Typography variant="body1">
          Add at least on filter and one fluorophore, and this tab will show a
          table of collection efficiency (sometimes called
          &quot;spillover&quot;) for each filter/fluorophore combination
        </Typography>
      </div>
    )
  }

  return (
    <div className="efficiency-table">
      <ErrorBoundary>
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
              onClick: () => setTransposed(transposed => !transposed),
            },
          ]}
        />
      </ErrorBoundary>
    </div>
  )
}

const OverlapToggle = ({ children, id, isActive }) => {
  const [active, setActive] = useState(isActive)
  const {
    data: { activeOverlaps },
  } = useQuery(GET_ACTIVE_OVERLAPS)

  useEffect(() => {
    setActive(activeOverlaps.includes(id))
  }, [activeOverlaps, id])

  const [setOverlaps] = useMutation(gql`
    mutation updateActiveOverlaps($add: [String], $remove: [String]) {
      updateActiveOverlaps(add: $add, remove: $remove) @client
    }
  `)

  const handleClick = event => {
    const elem = event.target.closest("button")
    const checked = !elem.checked
    const variables = {}
    variables[checked ? "add" : "remove"] = [id]
    setOverlaps({ variables })
    setActive(checked)
  }

  return (
    <Button
      size="small"
      onClick={handleClick}
      value={id}
      checked={active}
      variant={active ? "contained" : "outlined"}
      color={
        +children > 50 ? "primary" : +children > 5 ? "default" : "secondary"
      }
    >
      {children}
    </Button>
  )
}

export default EfficiencyTable

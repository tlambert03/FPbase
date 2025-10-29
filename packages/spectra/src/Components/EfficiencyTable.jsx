import { useApolloClient, useMutation, useQuery } from "@apollo/client"
import SaveAlt from "@mui/icons-material/SaveAlt"
import Shuffle from "@mui/icons-material/Shuffle"
import { Box, IconButton, Typography } from "@mui/material"
import Button from "@mui/material/Button"
import { makeStyles } from "@mui/styles"
import gql from "graphql-tag"
import { MaterialReactTable } from "material-react-table"
import React, { useEffect, useMemo, useState } from "react"
import { GET_ACTIVE_OVERLAPS } from "../client/queries"
import { trapz } from "../util"
import useSpectralData from "./useSpectraData"

class ErrorBoundary extends React.Component {
  static getDerivedStateFromError(_error) {
    return { hasError: true }
  }

  componentDidCatch(_error, _info) {}

  render() {
    const { children } = this.props
    return children
  }
}

const useStyles = makeStyles((_theme) => ({
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
    .map((arg) => arg.customId || arg.id)
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
  const [rows, setRows] = useState([])
  const classes = useStyles()
  const spectraData = useSpectralData()
  const client = useApolloClient()

  useEffect(() => {
    return () => {
      client.writeQuery({
        query: GET_ACTIVE_OVERLAPS,
        data: { activeOverlaps: [] },
      })
    }
  }, [client])

  // Generate table data
  useEffect(() => {
    async function updateTableData() {
      const filters = spectraData.filter(
        ({ category, subtype }) => category === "F" && subtype !== "BX"
      )
      const emSpectra = spectraData.filter(({ subtype }) => subtype === "EM")

      // untransposed columns represent different fluors
      let colItems = emSpectra
      let rowItems = filters
      const newRows = []
      if (transposed) {
        // columns represent different filters
        colItems = filters
        rowItems = emSpectra
      }

      rowItems.forEach((rowItem) => {
        const row = {
          field: rowItem.owner.name,
          _colItems: colItems,
          _rowItem: rowItem,
          _transposed: transposed,
        }
        colItems.forEach((colItem) => {
          const overlap = getOverlap(rowItem, colItem)
          const fluor = transposed ? rowItem : colItem
          row[colItem.owner.id] = ((100 * overlap.area) / fluor.area).toFixed(1)
          row[`${colItem.owner.id}_overlapID`] = overlap.id
        })
        newRows.push(row)
      })
      setRows(newRows)
    }

    updateTableData()
  }, [spectraData, transposed])

  // Generate columns using useMemo
  const columns = useMemo(() => {
    if (!rows.length) return []

    const filters = spectraData.filter(
      ({ category, subtype }) => category === "F" && subtype !== "BX"
    )
    const emSpectra = spectraData.filter(({ subtype }) => subtype === "EM")
    const colItems = transposed ? filters : emSpectra

    const cols = [
      {
        accessorKey: "field",
        header: transposed ? "Fluorophore" : "Filter",
        size: 150,
      },
    ]

    colItems.forEach(({ owner }) => {
      cols.push({
        accessorKey: owner.id,
        header: owner.name,
        size: 100,
        Cell: ({ row }) => {
          const overlapID = row.original[`${owner.id}_overlapID`]
          return <OverlapToggle id={overlapID}>{row.original[owner.id]}</OverlapToggle>
        },
      })
    })

    return cols
  }, [spectraData, transposed, rows.length])

  if (rows.length < 1 || columns.length < 2) {
    return (
      <div className={classes.description}>
        <Typography variant="h6">Efficiency Table</Typography>
        <Typography variant="body1">
          Add at least on filter and one fluorophore, and this tab will show a table of collection
          efficiency (sometimes called &quot;spillover&quot;) for each filter/fluorophore
          combination
        </Typography>
      </div>
    )
  }

  return (
    <Box className="efficiency-table">
      <ErrorBoundary>
        <MaterialReactTable
          columns={columns}
          data={rows}
          enablePagination={false}
          enableColumnActions={false}
          enableTopToolbar={true}
          enableBottomToolbar={false}
          enableSorting={true}
          enableFilters={false}
          enableGlobalFilter={false}
          muiTableProps={{
            sx: {
              tableLayout: "auto",
            },
          }}
          muiTableBodyCellProps={{
            sx: {
              fontSize: "1rem",
            },
          }}
          renderTopToolbarCustomActions={({ table }) => (
            <Box sx={{ display: "flex", gap: "1rem", alignItems: "center", p: 1 }}>
              <Typography variant="h6">Collection Efficiency (%)</Typography>
              <IconButton onClick={() => setTransposed((prev) => !prev)} title="Transpose">
                <Shuffle />
              </IconButton>
              <IconButton
                onClick={() => {
                  const csvContent = [
                    // Header row
                    columns
                      .map((col) => col.header)
                      .join(","),
                    // Data rows
                    ...rows.map((row) =>
                      columns.map((col) => row[col.accessorKey] || "").join(",")
                    ),
                  ].join("\n")

                  const blob = new Blob([csvContent], { type: "text/csv" })
                  const url = window.URL.createObjectURL(blob)
                  const a = document.createElement("a")
                  a.href = url
                  a.download = "efficiency-table.csv"
                  a.click()
                  window.URL.revokeObjectURL(url)
                }}
                title="Export to CSV"
              >
                <SaveAlt />
              </IconButton>
            </Box>
          )}
        />
      </ErrorBoundary>
    </Box>
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

  const handleClick = () => {
    const checked = !active
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
      color={+children > 50 ? "primary" : +children > 5 ? "inherit" : "secondary"}
    >
      {children}
    </Button>
  )
}

export default EfficiencyTable

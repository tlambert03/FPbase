import SaveAlt from "@mui/icons-material/SaveAlt"
import Shuffle from "@mui/icons-material/Shuffle"
import { Box, IconButton, Typography } from "@mui/material"
import Button from "@mui/material/Button"
import { makeStyles } from "@mui/styles"
import { MaterialReactTable } from "material-react-table"
import React, { useEffect, useMemo, useState } from "react"
import useSpectralData from "../hooks/useSpectraData"
import { useSpectraStore } from "../store/spectraStore"
import { computeOverlap } from "../utils/spectraUtils"

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

/**
 * Get or compute overlap spectrum with caching
 * Uses shared computeOverlap utility and spectraStore cache
 */
function getOverlap(store, ...spectra) {
  // Generate consistent ID for caching
  const sorted = [...spectra].sort((a, b) => {
    const aId = a.customId || a.id
    const bId = b.customId || b.id
    return String(aId).localeCompare(String(bId))
  })
  const idString = sorted.map((s) => s.customId || s.id).join("_")

  // Check cache first
  const cached = store.overlapCache[idString]
  if (cached) return cached

  // Compute and cache
  const overlap = computeOverlap(...spectra)
  store.setOverlapCache(idString, overlap)
  return overlap
}

const EfficiencyTable = ({ initialTranspose }) => {
  const [transposed, setTransposed] = useState(initialTranspose)
  const [rows, setRows] = useState([])
  const classes = useStyles()
  const spectraData = useSpectralData()
  const setActiveOverlaps = useSpectraStore((state) => state.setActiveOverlaps)
  const spectraStore = useSpectraStore()

  useEffect(() => {
    return () => {
      setActiveOverlaps([])
    }
  }, [setActiveOverlaps])

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
          const overlap = getOverlap(spectraStore, rowItem, colItem)
          const fluor = transposed ? rowItem : colItem
          row[colItem.owner.id] = ((100 * overlap.area) / fluor.area).toFixed(1)
          row[`${colItem.owner.id}_overlapID`] = overlap.id
        })
        newRows.push(row)
      })
      setRows(newRows)
    }

    updateTableData()
  }, [spectraData, transposed, spectraStore])

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
          renderTopToolbarCustomActions={() => (
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
  const activeOverlaps = useSpectraStore((state) => state.activeOverlaps)
  const updateActiveOverlaps = useSpectraStore((state) => state.updateActiveOverlaps)

  useEffect(() => {
    setActive(activeOverlaps.includes(id))
  }, [activeOverlaps, id])

  const handleClick = () => {
    const checked = !active
    if (checked) {
      updateActiveOverlaps([id])
    } else {
      updateActiveOverlaps(undefined, [id])
    }
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

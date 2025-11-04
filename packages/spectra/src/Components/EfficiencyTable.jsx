import SaveAlt from "@mui/icons-material/SaveAlt"
import Shuffle from "@mui/icons-material/Shuffle"
import {
  Box,
  IconButton,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TableSortLabel,
  Toolbar,
  Typography,
} from "@mui/material"
import Button from "@mui/material/Button"
import { makeStyles } from "@mui/styles"
import React, { useEffect, useMemo, useState } from "react"
import useSpectralData from "../hooks/useSpectraData"
import { useSpectraStore } from "../store/spectraStore"
import { computeOverlap, sortSpectraById } from "../utils/spectraUtils"

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
  toolbar: {
    display: "flex",
    gap: "1rem",
    alignItems: "center",
    padding: "8px 16px",
  },
}))

/**
 * Get or compute overlap spectrum with caching
 * Uses shared computeOverlap utility and spectraStore cache
 */
function getOverlap(store, ...spectra) {
  // Generate consistent ID for caching
  const sorted = sortSpectraById(spectra)
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
  const [orderBy, setOrderBy] = useState("field")
  const [order, setOrder] = useState("asc")
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
      })
    })

    return cols
  }, [spectraData, transposed, rows.length])

  // Sorting logic
  const handleSort = (columnId) => {
    const isAsc = orderBy === columnId && order === "asc"
    setOrder(isAsc ? "desc" : "asc")
    setOrderBy(columnId)
  }

  const sortedRows = useMemo(() => {
    if (!orderBy) return rows

    return [...rows].sort((a, b) => {
      const aValue = a[orderBy]
      const bValue = b[orderBy]

      // Handle numeric values
      const aNum = Number.parseFloat(aValue)
      const bNum = Number.parseFloat(bValue)
      if (!Number.isNaN(aNum) && !Number.isNaN(bNum)) {
        return order === "asc" ? aNum - bNum : bNum - aNum
      }

      // Handle string values
      const aStr = String(aValue || "")
      const bStr = String(bValue || "")
      const comparison = aStr.localeCompare(bStr)
      return order === "asc" ? comparison : -comparison
    })
  }, [rows, orderBy, order])

  // CSV Export
  const handleExportCSV = () => {
    const csvContent = [
      // Header row
      columns
        .map((col) => col.header)
        .join(","),
      // Data rows
      ...sortedRows.map((row) => columns.map((col) => row[col.accessorKey] || "").join(",")),
    ].join("\n")

    const blob = new Blob([csvContent], { type: "text/csv" })
    const url = window.URL.createObjectURL(blob)
    const a = document.createElement("a")
    a.href = url
    a.download = "efficiency-table.csv"
    a.click()
    window.URL.revokeObjectURL(url)
  }

  if (rows.length < 1 || columns.length < 2) {
    return (
      <div className={classes.description}>
        <Typography variant="h6">Efficiency Table</Typography>
        <Typography variant="body1">
          Add at least one filter and one fluorophore, and this tab will show a table of collection
          efficiency (sometimes called "spillover") for each filter/fluorophore combination.
        </Typography>
      </div>
    )
  }

  return (
    <Box className="efficiency-table">
      <ErrorBoundary>
        <Paper>
          <Toolbar className={classes.toolbar}>
            <Typography variant="h6">Collection Efficiency (%)</Typography>
            <IconButton onClick={() => setTransposed((prev) => !prev)} title="Transpose">
              <Shuffle />
            </IconButton>
            <IconButton onClick={handleExportCSV} title="Export to CSV">
              <SaveAlt />
            </IconButton>
          </Toolbar>

          <TableContainer>
            <Table className={classes.table} sx={{ tableLayout: "auto" }} size="small">
              <TableHead>
                <TableRow>
                  {columns.map((column) => (
                    <TableCell key={column.accessorKey} sx={{ width: column.size }}>
                      <TableSortLabel
                        active={orderBy === column.accessorKey}
                        direction={orderBy === column.accessorKey ? order : "asc"}
                        onClick={() => handleSort(column.accessorKey)}
                      >
                        {column.header}
                      </TableSortLabel>
                    </TableCell>
                  ))}
                </TableRow>
              </TableHead>
              <TableBody>
                {sortedRows.map((row) => (
                  <TableRow key={row.field} hover>
                    {columns.map((column) => (
                      <TableCell key={column.accessorKey} sx={{ fontSize: "1rem" }}>
                        {column.accessorKey === "field" ? (
                          row[column.accessorKey]
                        ) : (
                          <OverlapToggle id={row[`${column.accessorKey}_overlapID`]}>
                            {row[column.accessorKey]}
                          </OverlapToggle>
                        )}
                      </TableCell>
                    ))}
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          </TableContainer>
        </Paper>
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

  const value = Number(children)
  const color = value > 50 ? "primary" : value > 5 ? "inherit" : "secondary"

  return (
    <Button
      size="small"
      onClick={handleClick}
      value={id}
      checked={active}
      variant={active ? "contained" : "outlined"}
      color={color}
    >
      {children}
    </Button>
  )
}

export default EfficiencyTable

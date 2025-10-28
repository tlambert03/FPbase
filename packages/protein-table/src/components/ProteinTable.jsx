import AddCircleIcon from "@mui/icons-material/AddCircle"
import ArrowDownwardIcon from "@mui/icons-material/ArrowDownward"
import ArrowUpwardIcon from "@mui/icons-material/ArrowUpward"
import ChevronLeftIcon from "@mui/icons-material/ChevronLeft"
import ChevronRightIcon from "@mui/icons-material/ChevronRight"
import HelpOutlineIcon from "@mui/icons-material/HelpOutline"
import {
  alpha,
  Box,
  IconButton,
  MenuItem,
  Paper,
  Select,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Typography,
} from "@mui/material"
import {
  flexRender,
  getCoreRowModel,
  getPaginationRowModel,
  getSortedRowModel,
  useReactTable,
} from "@tanstack/react-table"
import { useMemo, useState } from "react"

/**
 * Row vertical padding - controls table row height
 * 0.5 = compact, 0.75 = comfortable, 1 = spacious
 */
const ROW_VERTICAL_PADDING = 0.7

/**
 * Wavelength thresholds for determining text color on colored backgrounds
 * Wavelengths between these values appear yellowish and need dark text
 */
const WAVELENGTH_DARK_TEXT_MIN = 477
const WAVELENGTH_DARK_TEXT_MAX = 590

/**
 * Opacity for colored spectral cells (ex_max, em_max)
 */
const SPECTRAL_CELL_OPACITY = 0.7

/**
 * Alpha values for hover effects
 */
const TABLE_HOVER_ALPHA = 0.08
const HEADER_HOVER_ALPHA = 0.1

/**
 * Glossary URL mapping for help links
 */
const GLOSSARY_BASE = "https://help.fpbase.org/glossary"
const GLOSSARY_LINKS = {
  ex_max: `${GLOSSARY_BASE}#ex-max`,
  em_max: `${GLOSSARY_BASE}#em-max`,
  stokes: `${GLOSSARY_BASE}#stokes-shift`,
  ext_coeff: `${GLOSSARY_BASE}#extinction-coefficient`,
  qy: `${GLOSSARY_BASE}#quantum-yield`,
  brightness: `${GLOSSARY_BASE}#brightness`,
  pka: `${GLOSSARY_BASE}#pka`,
  agg: `${GLOSSARY_BASE}#oligomerization`,
  maturation: `${GLOSSARY_BASE}#maturation`,
  lifetime: `${GLOSSARY_BASE}#lifetime`,
  weight: `${GLOSSARY_BASE}#molecular-weight`,
}

/**
 * Helper component for column header with optional help link
 */
function ColumnHeader({ children, glossaryKey }) {
  if (!glossaryKey || !GLOSSARY_LINKS[glossaryKey]) {
    return <>{children}</>
  }

  return (
    <Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
      {children}
      <a
        href={GLOSSARY_LINKS[glossaryKey]}
        target="_blank"
        rel="noopener noreferrer"
        onClick={(e) => e.stopPropagation()}
        style={{
          color: "inherit",
          fontSize: "0.75rem",
          opacity: 0.5,
          display: "flex",
          alignItems: "center",
          textDecoration: "none",
        }}
        onMouseEnter={(e) => {
          e.currentTarget.style.opacity = 0.8
          e.currentTarget.style.cursor = "help"
        }}
        onMouseLeave={(e) => {
          e.currentTarget.style.opacity = 0.5
        }}
      >
        <HelpOutlineIcon fontSize="inherit" />
      </a>
    </Box>
  )
}

/**
 * Format number with thousand separators
 */
function formatNumber(num) {
  if (num == null) return ""
  return num.toLocaleString()
}

/**
 * Reusable component for spectral wavelength cells with colored backgrounds
 */
function ColoredCell({ value, backgroundColor }) {
  if (!value) return ""

  const textColor =
    value > WAVELENGTH_DARK_TEXT_MIN && value < WAVELENGTH_DARK_TEXT_MAX ? "#000" : "#eee"

  return (
    <Box
      sx={{
        backgroundColor: backgroundColor || "#fff",
        color: textColor,
        opacity: SPECTRAL_CELL_OPACITY,
        display: "flex",
        alignItems: "center",
        justifyContent: "center",
        width: "100%",
        py: ROW_VERTICAL_PADDING,
      }}
    >
      {value}
    </Box>
  )
}

/**
 * Create column definitions for the protein table
 * Uses accessor functions to properly access nested data for sorting
 */
function createColumns() {
  return [
    {
      id: "name",
      accessorFn: (row) => row.protein.name,
      header: "Name",
      cell: ({ row }) => {
        const protein = row.original.protein
        const state = row.original.state
        return (
          <Box>
            <Box
              component="a"
              href={protein.url}
              sx={{
                textDecoration: "none",
                color: "primary.main",
                fontWeight: 500,
              }}
            >
              {protein.name}
            </Box>
            {state && state.name !== "default" && (
              <Box
                component="span"
                sx={{
                  ml: 0.5,
                  fontSize: "0.85em",
                  color: "text.secondary",
                }}
              >
                ({state.name})
              </Box>
            )}
          </Box>
        )
      },
      size: 220,
    },
    {
      id: "ex_max",
      accessorFn: (row) => row.state?.ex_max ?? null,
      header: () => <ColumnHeader glossaryKey="ex_max">λₑₓ</ColumnHeader>,
      cell: ({ row }) => {
        const state = row.original.state
        return <ColoredCell value={state?.ex_max} backgroundColor={state?.exhex} />
      },
      size: 50,
      meta: { align: "center", noPadding: true },
    },
    {
      id: "em_max",
      accessorFn: (row) => row.state?.em_max ?? null,
      header: () => <ColumnHeader glossaryKey="em_max">λₑₘ</ColumnHeader>,
      cell: ({ row }) => {
        const state = row.original.state
        return <ColoredCell value={state?.em_max} backgroundColor={state?.emhex} />
      },
      size: 50,
      meta: { align: "center", noPadding: true },
    },
    {
      id: "stokes",
      accessorFn: (row) => row.state?.stokes ?? null,
      header: () => <ColumnHeader glossaryKey="stokes">Stokes</ColumnHeader>,
      cell: ({ row }) => row.original.state?.stokes || "",
      size: 40,
      meta: { align: "center" },
    },
    {
      id: "ext_coeff",
      accessorFn: (row) => row.state?.ext_coeff ?? null,
      header: () => <ColumnHeader glossaryKey="ext_coeff">EC</ColumnHeader>,
      cell: ({ row }) => formatNumber(row.original.state?.ext_coeff),
      size: 60,
      meta: { align: "right" },
    },
    {
      id: "qy",
      accessorFn: (row) => row.state?.qy ?? null,
      header: () => <ColumnHeader glossaryKey="qy">QY</ColumnHeader>,
      cell: ({ row }) => row.original.state?.qy || "",
      size: 60,
      meta: { align: "center" },
    },
    {
      id: "brightness",
      accessorFn: (row) => row.state?.brightness ?? null,
      header: () => <ColumnHeader glossaryKey="brightness">Brightness</ColumnHeader>,
      cell: ({ row }) => row.original.state?.brightness || "",
      size: 60,
      meta: { align: "center" },
    },
    {
      id: "pka",
      accessorFn: (row) => row.state?.pka ?? null,
      header: () => <ColumnHeader glossaryKey="pka">pKa</ColumnHeader>,
      cell: ({ row }) => row.original.state?.pka || "",
      size: 50,
      meta: { align: "center" },
    },
    {
      id: "agg",
      accessorFn: (row) => row.protein.agg ?? "",
      header: () => <ColumnHeader glossaryKey="agg">Agg</ColumnHeader>,
      cell: ({ row }) => row.original.protein.agg || "",
      size: 50,
      meta: { align: "center" },
    },
    {
      id: "maturation",
      accessorFn: (row) => row.state?.maturation ?? null,
      header: () => <ColumnHeader glossaryKey="maturation">Mat</ColumnHeader>,
      cell: ({ row }) => row.original.state?.maturation || "",
      size: 50,
      meta: { align: "center" },
    },
    {
      id: "lifetime",
      accessorFn: (row) => row.state?.lifetime ?? null,
      header: () => <ColumnHeader glossaryKey="lifetime">τ</ColumnHeader>,
      cell: ({ row }) => row.original.state?.lifetime || "",
      size: 50,
      meta: { align: "center" },
    },
    {
      id: "weight",
      accessorFn: (row) => row.protein.weight ?? null,
      header: () => <ColumnHeader glossaryKey="weight">kDa</ColumnHeader>,
      cell: ({ row }) => row.original.protein.weight || "",
      size: 70,
      meta: { align: "center" },
    },
    {
      id: "year",
      accessorFn: (row) => row.protein.year ?? null,
      header: "Year",
      cell: ({ row }) => row.original.protein.year || "",
      size: 60,
      meta: { align: "center" },
    },
    {
      id: "compare",
      header: () => <Box sx={{ display: { xs: "none", lg: "block" } }}>Compare</Box>,
      cell: ({ row }) => (
        <IconButton
          size="small"
          className="comparison-btn"
          data-flash="1"
          data-action-url="/ajax/comparison/"
          data-object={row.original.protein.slug}
          data-op="add"
          sx={{
            color: "primary.main",
            p: 0.5,
            "&:hover": {
              backgroundColor: (theme) => alpha(theme.palette.primary.main, TABLE_HOVER_ALPHA),
            },
          }}
        >
          <AddCircleIcon sx={{ fontSize: "1rem", color: "gray" }} />
        </IconButton>
      ),
      size: 60,
      meta: { align: "center", noPadding: true },
      enableSorting: false,
    },
  ]
}

/**
 * Transform protein data to flat rows (one row per non-dark state)
 */
function transformProteinsToRows(proteins) {
  const rows = []

  proteins.forEach((protein) => {
    const states = protein.states?.filter((s) => !s.is_dark) || []

    if (states.length === 0) {
      // Protein with no states
      rows.push({ protein, state: null })
    } else {
      // Create one row per state
      states.forEach((state) => {
        rows.push({ protein, state })
      })
    }
  })

  return rows
}

/**
 * Main ProteinTable component using TanStack Table
 */
export default function ProteinTable({ proteins, filters, totalCount }) {
  const [sorting, setSorting] = useState([{ id: "brightness", desc: true }])
  const [pagination, setPagination] = useState({
    pageIndex: 0,
    pageSize: 25,
  })

  // Transform proteins to flat rows
  const allData = useMemo(() => transformProteinsToRows(proteins), [proteins])

  // Apply filters from external filter controls
  const filteredData = useMemo(() => {
    if (!filters.switchType && !filters.aggType && !filters.search) return allData

    return allData.filter((row) => {
      const protein = row.protein

      // Search filter (case-insensitive protein name)
      if (filters.search) {
        const searchLower = filters.search.toLowerCase()
        const nameMatch = protein.name.toLowerCase().includes(searchLower)
        if (!nameMatch) return false
      }

      // Switch type filter
      if (filters.switchType && protein.switch_type !== filters.switchType) {
        return false
      }

      // Aggregation filter
      if (filters.aggType && protein.agg !== filters.aggType) {
        return false
      }

      return true
    })
  }, [allData, filters])

  const columns = useMemo(() => createColumns(), [])

  const table = useReactTable({
    data: filteredData,
    columns,
    state: {
      sorting,
      pagination,
    },
    onSortingChange: setSorting,
    onPaginationChange: setPagination,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getPaginationRowModel: getPaginationRowModel(),
    manualPagination: false, // Client-side pagination
  })

  const handleChangePage = (newPage) => {
    setPagination((old) => ({ ...old, pageIndex: newPage }))
  }

  const handleChangeRowsPerPage = (event) => {
    setPagination({
      pageIndex: 0,
      pageSize: parseInt(event.target.value, 10),
    })
  }

  const canPreviousPage = pagination.pageIndex > 0
  const canNextPage =
    pagination.pageIndex < Math.ceil(filteredData.length / pagination.pageSize) - 1
  const startRow = pagination.pageIndex * pagination.pageSize + 1
  const endRow = Math.min((pagination.pageIndex + 1) * pagination.pageSize, filteredData.length)

  const isFiltered = filters.switchType || filters.aggType || filters.search
  const displayText = isFiltered
    ? `Showing ${
        table.getRowModel().rows.length
      } of ${filteredData.length} entries (filtered from ${totalCount} total)`
    : `Showing ${table.getRowModel().rows.length} of ${filteredData.length} entries`

  return (
    <Box>
      <TableContainer component={Paper} sx={{ boxShadow: 2 }}>
        <Table size="small" sx={{ minWidth: 960 }}>
          <TableHead>
            {table.getHeaderGroups().map((headerGroup) => (
              <TableRow key={headerGroup.id}>
                {headerGroup.headers.map((header) => (
                  <TableCell
                    key={header.id}
                    align={header.column.columnDef.meta?.align || "left"}
                    sx={{
                      fontWeight: 600,
                      fontSize: "0.75rem",
                      backgroundColor: "grey.100",
                      cursor: header.column.getCanSort() ? "pointer" : "default",
                      userSelect: "none",
                      width: header.column.columnDef.size,
                      py: 0.5,
                      px: 1,
                      "&:hover": header.column.getCanSort()
                        ? {
                            backgroundColor: alpha("#000", HEADER_HOVER_ALPHA),
                          }
                        : {},
                    }}
                    onClick={header.column.getToggleSortingHandler()}
                  >
                    <Box
                      sx={{
                        display: "flex",
                        alignItems: "center",
                        justifyContent:
                          header.column.columnDef.meta?.align === "right"
                            ? "flex-end"
                            : header.column.columnDef.meta?.align === "center"
                              ? "center"
                              : "flex-start",
                        gap: 0.5,
                      }}
                    >
                      {flexRender(header.column.columnDef.header, header.getContext())}
                      {header.column.getIsSorted() && (
                        <Box
                          sx={{
                            display: "flex",
                            alignItems: "center",
                            fontSize: "0.75rem",
                          }}
                        >
                          {header.column.getIsSorted() === "desc" ? (
                            <ArrowDownwardIcon fontSize="inherit" />
                          ) : (
                            <ArrowUpwardIcon fontSize="inherit" />
                          )}
                        </Box>
                      )}
                    </Box>
                  </TableCell>
                ))}
              </TableRow>
            ))}
          </TableHead>
          <TableBody>
            {table.getRowModel().rows.map((row) => (
              <TableRow
                key={row.id}
                sx={{
                  "&:nth-of-type(odd)": {
                    backgroundColor: "action.hover",
                  },
                  "&:hover": {
                    backgroundColor: (theme) =>
                      alpha(theme.palette.primary.main, TABLE_HOVER_ALPHA),
                  },
                }}
              >
                {row.getVisibleCells().map((cell) => {
                  const isColorCell = cell.column.columnDef.meta?.noPadding
                  return (
                    <TableCell
                      key={cell.id}
                      align={cell.column.columnDef.meta?.align || "left"}
                      sx={{
                        py: isColorCell ? 0 : ROW_VERTICAL_PADDING,
                        px: isColorCell ? 0 : 0.5,
                        fontSize: "0.75rem",
                      }}
                    >
                      {flexRender(cell.column.columnDef.cell, cell.getContext())}
                    </TableCell>
                  )
                })}
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </TableContainer>

      <Paper
        sx={{
          mt: 2,
          px: 2,
          py: 1.5,
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
        }}
      >
        <Typography variant="body2" color="text.secondary">
          {displayText}
        </Typography>
        <Box sx={{ display: "flex", alignItems: "center", gap: 2 }}>
          <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
            <Typography variant="body2" color="text.secondary">
              Rows per page:
            </Typography>
            <Select
              value={pagination.pageSize}
              onChange={handleChangeRowsPerPage}
              size="small"
              sx={{
                fontSize: "0.875rem",
                "& .MuiSelect-select": {
                  py: 0.5,
                  pr: 3,
                },
              }}
            >
              <MenuItem value={10}>10</MenuItem>
              <MenuItem value={25}>25</MenuItem>
              <MenuItem value={50}>50</MenuItem>
              <MenuItem value={100}>100</MenuItem>
            </Select>
          </Box>
          <Typography variant="body2" color="text.secondary">
            {startRow}–{endRow} of {filteredData.length}
          </Typography>
          <Box sx={{ display: "flex", alignItems: "center" }}>
            <IconButton
              onClick={() => handleChangePage(pagination.pageIndex - 1)}
              disabled={!canPreviousPage}
              size="small"
            >
              <ChevronLeftIcon />
            </IconButton>
            <IconButton
              onClick={() => handleChangePage(pagination.pageIndex + 1)}
              disabled={!canNextPage}
              size="small"
            >
              <ChevronRightIcon />
            </IconButton>
          </Box>
        </Box>
      </Paper>
    </Box>
  )
}

import DownloadIcon from "@mui/icons-material/Download"
import SearchIcon from "@mui/icons-material/Search"
import {
  Box,
  Button,
  FormControl,
  InputLabel,
  MenuItem,
  Paper,
  Select,
  TextField,
} from "@mui/material"
import { exportToCSV, prepareExportData } from "../utils/export"

const SWITCH_TYPES = [
  { value: "", label: "All" },
  { value: "b", label: "Basic" },
  { value: "pa", label: "Photoactivatable" },
  { value: "ps", label: "Photoswitchable" },
  { value: "pc", label: "Photoconvertible" },
  { value: "mp", label: "Multi-photochromic" },
  { value: "t", label: "Timer" },
  { value: "o", label: "Other" },
]

const AGG_TYPES = [
  { value: "", label: "All" },
  { value: "m", label: "Monomer" },
  { value: "d", label: "Dimer" },
  { value: "td", label: "Tandem dimer" },
  { value: "wd", label: "Weak dimer" },
  { value: "t", label: "Tetramer" },
]

/**
 * Table control bar with filters and export buttons
 */
export default function TableControls({ proteins, filters, onFilterChange }) {
  const handleExportCSV = () => {
    const data = prepareExportData(proteins)
    exportToCSV(data, "fpbase_proteins.csv")
  }

  return (
    <Paper
      elevation={1}
      sx={{
        p: 2,
        mb: 2,
        display: "flex",
        flexWrap: "wrap",
        gap: 2,
        alignItems: "center",
      }}
    >
      <FormControl size="small" sx={{ minWidth: 200 }}>
        <InputLabel>Switch Type</InputLabel>
        <Select
          value={filters.switchType}
          label="Switch Type"
          onChange={(e) => onFilterChange({ ...filters, switchType: e.target.value })}
        >
          {SWITCH_TYPES.map((type) => (
            <MenuItem key={type.value} value={type.value}>
              {type.label}
            </MenuItem>
          ))}
        </Select>
      </FormControl>

      <FormControl size="small" sx={{ minWidth: 200 }}>
        <InputLabel>Oligomerization</InputLabel>
        <Select
          value={filters.aggType}
          label="Oligomerization"
          onChange={(e) => onFilterChange({ ...filters, aggType: e.target.value })}
        >
          {AGG_TYPES.map((type) => (
            <MenuItem key={type.value} value={type.value}>
              {type.label}
            </MenuItem>
          ))}
        </Select>
      </FormControl>

      <TextField
        size="small"
        placeholder="Search proteins..."
        value={filters.search || ""}
        onChange={(e) => onFilterChange({ ...filters, search: e.target.value })}
        InputProps={{
          startAdornment: <SearchIcon sx={{ mr: 1, color: "text.secondary" }} />,
        }}
        sx={{ minWidth: 250 }}
      />

      <Box sx={{ flexGrow: 1 }} />

      <Button
        variant="outlined"
        size="small"
        startIcon={<DownloadIcon />}
        onClick={handleExportCSV}
      >
        Export CSV
      </Button>
    </Paper>
  )
}

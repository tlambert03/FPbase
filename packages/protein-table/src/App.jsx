import { useState, useEffect } from "react"
import {
  Box,
  Container,
  Typography,
  CircularProgress,
  Alert,
  ThemeProvider,
  createTheme,
} from "@mui/material"
import ProteinTable from "./components/ProteinTable"
import TableControls from "./components/TableControls"

const theme = createTheme({
  palette: {
    primary: {
      main: "#0066cc",
    },
  },
  typography: {
    fontFamily: [
      "-apple-system",
      "BlinkMacSystemFont",
      '"Segoe UI"',
      "Roboto",
      '"Helvetica Neue"',
      "Arial",
      "sans-serif",
    ].join(","),
  },
})

/**
 * Main App component for the protein table viewer
 */
export default function App() {
  const [proteins, setProteins] = useState([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [filters, setFilters] = useState({
    search: "",
    switchType: "",
    aggType: "",
  })

  useEffect(() => {
    async function fetchProteins() {
      try {
        const response = await fetch("/api/proteins/table-data/")
        if (!response.ok) {
          throw new Error(`HTTP error! status: ${response.status}`)
        }
        const data = await response.json()
        setProteins(data)
        setLoading(false)
      } catch (err) {
        console.error("Failed to fetch proteins:", err)
        setError(err.message)
        setLoading(false)
      }
    }

    fetchProteins()
  }, [])

  if (loading) {
    return (
      <ThemeProvider theme={theme}>
        <Container maxWidth="xl" sx={{ py: 4 }}>
          <Box
            sx={{
              display: "flex",
              justifyContent: "center",
              alignItems: "center",
              minHeight: 400,
            }}
          >
            <CircularProgress />
          </Box>
        </Container>
      </ThemeProvider>
    )
  }

  if (error) {
    return (
      <ThemeProvider theme={theme}>
        <Container maxWidth="xl" sx={{ py: 4 }}>
          <Alert severity="error">
            Failed to load protein data: {error}
          </Alert>
        </Container>
      </ThemeProvider>
    )
  }

  return (
    <ThemeProvider theme={theme}>
      <Container maxWidth="xl" sx={{ py: 3 }}>
        <Typography
          variant="h4"
          component="h2"
          align="center"
          gutterBottom
          sx={{ fontWeight: 600, mb: 2 }}
        >
          Fluorescent Protein Table
        </Typography>

        <Typography
          variant="body2"
          color="text.secondary"
          align="center"
          sx={{ mb: 3, maxWidth: 800, mx: "auto" }}
        >
          A comprehensive table of all proteins in the database. Click column
          headers to sort (shift-click for multi-sort). Use the filters below
          to narrow your search.
        </Typography>

        <TableControls
          proteins={proteins}
          filters={filters}
          onFilterChange={setFilters}
        />

        <ProteinTable
          proteins={proteins}
          filters={filters}
          totalCount={proteins.length}
        />
      </Container>
    </ThemeProvider>
  )
}

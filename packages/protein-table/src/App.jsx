import { Alert, Box, CircularProgress, Container, ThemeProvider } from '@mui/material'
import { QueryClient, QueryClientProvider, useQuery } from '@tanstack/react-query'
import { useState } from 'react'
import ProteinTable from './components/ProteinTable'
import TableControls from './components/TableControls'
import theme from './theme'

/**
 * Configure React Query client with appropriate caching settings
 */
const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      staleTime: 5 * 60 * 1000, // 5 minutes
      gcTime: 30 * 60 * 1000, // 30 minutes (previously cacheTime)
      refetchOnWindowFocus: false,
      retry: 1,
    },
  },
})

/**
 * Fetch proteins from the API
 */
async function fetchProteins() {
  const response = await fetch('/api/proteins/table-data/')
  if (!response.ok) {
    throw new Error(`HTTP error! status: ${response.status}`)
  }
  return response.json()
}

/**
 * Inner App component that uses React Query hooks
 */
function AppContent() {
  const [filters, setFilters] = useState({
    search: '',
    switchType: '',
    aggType: '',
  })

  const {
    data: proteins = [],
    isLoading,
    error,
  } = useQuery({
    queryKey: ['proteins', 'table-data'],
    queryFn: fetchProteins,
    onError: (error) => {
      if (process.env.NODE_ENV === 'development') {
        console.error('Failed to fetch proteins:', error)
      }
    },
  })

  return (
    <ThemeProvider theme={theme}>
      <Container maxWidth="xl" sx={{ py: isLoading || error ? 4 : 3 }}>
        {isLoading && (
          <Box
            sx={{
              display: 'flex',
              justifyContent: 'center',
              alignItems: 'center',
              minHeight: 400,
            }}
          >
            <CircularProgress />
          </Box>
        )}

        {error && <Alert severity="error">Failed to load protein data: {error.message}</Alert>}

        {!isLoading && !error && (
          <>
            <TableControls proteins={proteins} filters={filters} onFilterChange={setFilters} />

            <ProteinTable proteins={proteins} filters={filters} totalCount={proteins.length} />
          </>
        )}
      </Container>
    </ThemeProvider>
  )
}

/**
 * Main App component wrapped with React Query provider
 */
export default function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <AppContent />
    </QueryClientProvider>
  )
}

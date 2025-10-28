import { createTheme } from '@mui/material'

/**
 * Material-UI theme configuration for the protein table
 */
const theme = createTheme({
  palette: {
    primary: {
      main: '#0066cc',
    },
  },
  typography: {
    fontFamily: [
      '-apple-system',
      'BlinkMacSystemFont',
      '"Segoe UI"',
      'Roboto',
      '"Helvetica Neue"',
      'Arial',
      'sans-serif',
    ].join(','),
  },
})

export default theme

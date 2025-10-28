import { createTheme } from '@mui/material/styles'

const theme = createTheme({
  //   palette: {
  //     primary: purple,
  //     secondary: green
  //   },
  //   status: {
  //     danger: "orange"
  //   },
  breakpoints: {
    values: {
      xs: 0,
      sm: 768,
      md: 960,
      lg: 1280,
      sl: 1920,
    },
  },
})

export default theme

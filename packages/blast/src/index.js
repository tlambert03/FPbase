import React from "react"
import { ThemeProvider, StyledEngineProvider, createTheme } from "@mui/material/styles"
import App from "./App"

const theme = createTheme()

const AppWrapper = () => {
  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <App />
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

export default AppWrapper

import { createTheme, StyledEngineProvider, ThemeProvider } from "@mui/material/styles"
import React from "react"
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

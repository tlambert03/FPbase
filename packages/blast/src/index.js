/**
 * BLAST App Wrapper
 *
 * Uses the shared MUI theme configuration from design tokens
 * to ensure consistency with Bootstrap styling across the site.
 */

import { StyledEngineProvider, ThemeProvider, createTheme } from "@mui/material/styles"
import { getMuiThemeConfig } from "../../../frontend/src/theme/tokens.js"
import App from "./App"

const theme = createTheme(getMuiThemeConfig("light"))

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

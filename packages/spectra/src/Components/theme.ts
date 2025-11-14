/**
 * Spectra package theme
 *
 * Uses the shared MUI theme configuration from design tokens
 * to ensure consistency with Bootstrap styling across the site.
 */

import type { Theme } from "@mui/material/styles"
import { createTheme } from "@mui/material/styles"
import { getMuiThemeConfig } from "../../../../frontend/src/theme/tokens.js"

const theme: Theme = createTheme(getMuiThemeConfig("light"))

export default theme

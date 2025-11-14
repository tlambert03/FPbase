/**
 * Protein Table package theme
 *
 * Uses the shared MUI theme configuration from design tokens
 * to ensure consistency with Bootstrap styling across the site.
 *
 * Note: Previously this used #0066cc as the primary color, which was inconsistent
 * with the rest of the site. Now it uses the unified primary green (#3CA644).
 */

import { createTheme } from "@mui/material"
import { getMuiThemeConfig } from "../../../frontend/src/theme/tokens.js"

const theme = createTheme(getMuiThemeConfig("light"))

export default theme

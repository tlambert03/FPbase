/**
 * Protein Table package theme
 *
 * Creates MUI theme based on design tokens from frontend/src/theme/tokens.js
 * to ensure consistency with Bootstrap styling across the site.
 *
 * Note: Previously this used #0066cc as the primary color, which was inconsistent
 * with the rest of the site. Now it uses the unified primary green (#3CA644).
 */

import { createTheme } from "@mui/material"
import tokens from "../../../frontend/src/theme/tokens.js"

const theme = createTheme({
	palette: {
		mode: "light",
		primary: {
			main: tokens.colors.primary,
			dark: tokens.colors.primaryDark,
			light: tokens.colors.primaryLight,
		},
		secondary: {
			main: tokens.colors.secondary,
		},
		success: {
			main: tokens.colors.success,
		},
		error: {
			main: tokens.colors.danger,
		},
		warning: {
			main: tokens.colors.warning,
		},
		info: {
			main: tokens.colors.info,
		},
	},
	typography: {
		fontFamily: tokens.typography.fontFamily,
		h1: { fontWeight: tokens.typography.headingWeight },
		h2: { fontWeight: tokens.typography.headingWeight },
		h3: { fontWeight: tokens.typography.headingWeight },
		h4: { fontWeight: tokens.typography.headingWeight },
		h5: { fontWeight: tokens.typography.headingWeight },
		h6: { fontWeight: tokens.typography.headingWeight },
	},
	breakpoints: {
		values: {
			xs: 0,
			sm: 768,
			md: 960,
			lg: 1280,
			xl: 1920,
		},
	},
})

export default theme

/**
 * MUI Theme Generator - Creates Material-UI themes from design tokens
 *
 * This file creates a unified MUI theme based on the same design tokens
 * used for Bootstrap/SCSS styling. This ensures consistency between
 * the Django/SSR pages (Bootstrap) and React apps (MUI).
 *
 * Usage:
 *   import { createAppTheme } from 'path/to/muiTheme.js'
 *   const theme = createAppTheme('light')
 */

import { createTheme } from "@mui/material/styles"
import tokens from "./tokens.js"

/**
 * Creates a MUI theme from design tokens
 *
 * @param {('light'|'dark')} mode - The theme mode
 * @returns {import('@mui/material/styles').Theme} MUI theme object
 */
export function createAppTheme(mode = "light") {
	const modeTokens = tokens.modes[mode]

	return createTheme({
		palette: {
			mode,
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
			background: {
				default: modeTokens.background,
				paper: mode === "light" ? "#ffffff" : "#2a2a2a",
			},
			text: {
				primary: modeTokens.text,
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

		// Match Bootstrap breakpoints
		breakpoints: {
			values: {
				xs: 0,
				sm: 768, // Bootstrap sm breakpoint
				md: 960,
				lg: 1280,
				xl: 1920,
			},
		},

		// Component overrides for consistency with Bootstrap
		components: {
			MuiCssBaseline: {
				styleOverrides: `
          @import url('${tokens.typography.fontImportUrl}');

          :root {
            --navbar-height: ${tokens.spacing.navbarHeight};
            --color-primary: ${tokens.colors.primary};
            --color-footer-bg: ${tokens.colors.footerBg};
            --color-footer-text: ${tokens.colors.footerText};
          }
        `,
			},
		},
	})
}

/**
 * Default light theme
 */
export const lightTheme = createAppTheme("light")

/**
 * Default dark theme (for future use)
 */
export const darkTheme = createAppTheme("dark")

export default createAppTheme

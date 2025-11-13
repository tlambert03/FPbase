/**
 * Design Tokens - Single Source of Truth for FPbase Theming
 *
 * This file contains all design tokens (colors, typography, spacing, etc.)
 * used across both Bootstrap/SCSS and MUI/React parts of the application.
 *
 * These tokens are:
 * - Converted to SCSS variables (see generateScss.js)
 * - Used to generate MUI themes (see muiTheme.js)
 * - Can be converted to CSS custom properties for runtime theming
 */

export const tokens = {
	colors: {
		// Base color palette
		white: "#fff",
		black: "#000",

		// Grays
		gray100: "#f8f9fa",
		gray200: "#e9ecef",
		gray300: "#dee2e6",
		gray400: "#ced4da",
		gray500: "#adb5bd",
		gray600: "#868e96",
		gray700: "#495057",
		gray800: "#343a40",
		gray900: "#212529",

		// Brand colors (from Bootstrap variables)
		blue: "#3085d6",
		indigo: "#6610f2",
		purple: "#6f42c1",
		pink: "#e83e8c",
		red: "#dc3545",
		orange: "#fd7e14",
		yellow: "#ffc107",
		green: "#419d45",
		teal: "#20c997",
		cyan: "#17a2b8",

		// Theme colors
		// Note: Primary is saturate($green, 10%) in SCSS, which equals approximately #3CA644
		primary: "#3CA644",
		primaryDark: "#2E8935",
		primaryLight: "#5CB865",
		secondary: "#3085d6",
		success: "#419d45",
		danger: "#dc3545",
		warning: "#ffc107",
		info: "#17a2b8",
		light: "#f8f9fa",
		dark: "#343a40",

		// Site-specific colors
		footerBg: "#0C4B33",
		footerText: "#2B8C67",
		navbarBg: "#0D4B33", // Used in React components

		// Alert colors (from style.scss)
		alertDebugBg: "#fff",
		alertDebugBorder: "#d6e9c6",
		alertDebugText: "#000",
		alertErrorBg: "#f2dede",
		alertErrorBorder: "#eed3d7",
		alertErrorText: "#b94a48",

		// Other specific colors
		logoColor: "#efe",
	},

	typography: {
		// Font families
		fontFamily:
			"'Raleway', -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Roboto', 'Helvetica Neue', Arial, sans-serif",
		fontFamilyMonospace:
			"source-code-pro, Menlo, Monaco, Consolas, 'Courier New', monospace",

		// Font weights
		headingWeight: 100,
		bodyWeight: 400,

		// Font imports
		fontImportUrl: "https://fonts.googleapis.com/css?family=Raleway:100,500,600",
	},

	spacing: {
		navbarHeight: "3.5rem",
	},

	shadows: {
		enabled: true,
		dropdownShadow: "0 0.5rem 1rem rgba(0, 0, 0, 0.15)",
		navbarShadow: "0 0.1rem 0.2rem rgba(0, 0, 0, 0.2)",
	},

	options: {
		// Bootstrap options
		enableShadows: true,
		enableGradients: true,
		enableRounded: false, // Disabled by default in _options.scss
	},

	// Theme modes for future dark mode support
	modes: {
		light: {
			background: "#ffffff",
			text: "#212529",
			bodyBg: "#ffffff",
		},
		dark: {
			background: "#1a1a1a",
			text: "#f8f9fa",
			bodyBg: "#1a1a1a",
		},
	},
}

export default tokens

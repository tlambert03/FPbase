#!/usr/bin/env node

/**
 * Generate SCSS variables from design tokens
 *
 * This script reads tokens.js and generates _generated-variables.scss
 * which can be imported into the main SCSS files.
 *
 * Run this script before building the frontend to ensure SCSS has
 * the latest token values.
 */

import { writeFileSync } from "node:fs"
import { fileURLToPath } from "node:url"
import { dirname, join } from "node:path"
import tokens from "./tokens.js"

const __filename = fileURLToPath(import.meta.url)
const __dirname = dirname(__filename)

const outputPath = join(__dirname, "../css/_generated-variables.scss")

// Helper to create SCSS map syntax
function createScssMap(name, obj) {
	const entries = Object.entries(obj)
		.map(([key, value]) => `  "${key}": ${value}`)
		.join(",\n")
	return `$${name}: (\n${entries}\n);`
}

// Generate the SCSS file content
const scssContent = `// ============================================================================
// AUTO-GENERATED FILE - DO NOT EDIT DIRECTLY
// ============================================================================
// This file is generated from frontend/src/theme/tokens.js
// Run 'pnpm generate-scss' to regenerate
// ============================================================================

// Font imports
@import url('${tokens.typography.fontImportUrl}');

// ============================================================================
// Color System
// ============================================================================

// Base colors
$white:    ${tokens.colors.white};
$black:    ${tokens.colors.black};

// Grayscale
$gray-100: ${tokens.colors.gray100};
$gray-200: ${tokens.colors.gray200};
$gray-300: ${tokens.colors.gray300};
$gray-400: ${tokens.colors.gray400};
$gray-500: ${tokens.colors.gray500};
$gray-600: ${tokens.colors.gray600};
$gray-700: ${tokens.colors.gray700};
$gray-800: ${tokens.colors.gray800};
$gray-900: ${tokens.colors.gray900};

$grays: (
  "100": $gray-100,
  "200": $gray-200,
  "300": $gray-300,
  "400": $gray-400,
  "500": $gray-500,
  "600": $gray-600,
  "700": $gray-700,
  "800": $gray-800,
  "900": $gray-900
);

// Brand colors
$blue:    ${tokens.colors.blue};
$indigo:  ${tokens.colors.indigo};
$purple:  ${tokens.colors.purple};
$pink:    ${tokens.colors.pink};
$red:     ${tokens.colors.red};
$orange:  ${tokens.colors.orange};
$yellow:  ${tokens.colors.yellow};
$green:   ${tokens.colors.green};
$teal:    ${tokens.colors.teal};
$cyan:    ${tokens.colors.cyan};

$colors: (
  "blue":       $blue,
  "indigo":     $indigo,
  "purple":     $purple,
  "pink":       $pink,
  "red":        $red,
  "orange":     $orange,
  "yellow":     $yellow,
  "green":      $green,
  "teal":       $teal,
  "cyan":       $cyan,
  "white":      $white,
  "gray":       $gray-600,
  "gray-dark":  $gray-800
);

// Theme colors (Bootstrap semantic colors)
$primary:       ${tokens.colors.primary};
$secondary:     ${tokens.colors.secondary};
$success:       ${tokens.colors.success};
$info:          ${tokens.colors.info};
$warning:       ${tokens.colors.warning};
$danger:        ${tokens.colors.danger};
$light:         ${tokens.colors.light};
$dark:          ${tokens.colors.dark};

$theme-colors: (
  "primary":    $primary,
  "secondary":  $secondary,
  "success":    $success,
  "info":       $info,
  "warning":    $warning,
  "danger":     $danger,
  "light":      $light,
  "dark":       $dark
);

// ============================================================================
// Site-specific Colors
// ============================================================================

$footer-background: ${tokens.colors.footerBg};
$footer-color:      ${tokens.colors.footerText};
$navbar-bg:         ${tokens.colors.navbarBg};
$logo-color:        ${tokens.colors.logoColor};

// ============================================================================
// Typography
// ============================================================================

$headings-font-weight: ${tokens.typography.headingWeight};
$body-font-weight:     ${tokens.typography.bodyWeight};

// Note: Font family is loaded via @import url() above
// Bootstrap will use its default font stack unless overridden

// ============================================================================
// Spacing
// ============================================================================

$navbar-height: ${tokens.spacing.navbarHeight};

// ============================================================================
// Component Customization
// ============================================================================

// Dropdowns
$dropdown-bg: $white;
$dropdown-link-hover-bg: rgba($primary, 0.3);
$dropdown-box-shadow: ${tokens.shadows.dropdownShadow};

// Navbar
$navbar-dark-color:        rgba($white, 0.8);
$navbar-dark-hover-color:  rgba($white, 0.95);

// Tooltips
$tooltip-text: $white;
$tooltip-color: rgba($footer-background, 0.8);
`

try {
	writeFileSync(outputPath, scssContent, "utf-8")
	console.log("✓ Generated SCSS variables from tokens")
	console.log(`  Output: ${outputPath}`)
} catch (error) {
	console.error("✗ Failed to generate SCSS variables:", error)
	process.exit(1)
}

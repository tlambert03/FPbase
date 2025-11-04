/**
 * Type augmentation for @mui/styles to use MUI theme
 */
import type { Theme } from "@mui/material/styles"

declare module "@mui/styles" {
  // eslint-disable-next-line @typescript-eslint/no-empty-object-type
  interface DefaultTheme extends Theme {}
}

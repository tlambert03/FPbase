import react from "@vitejs/plugin-react"
import { defineConfig } from "vitest/config"

export default defineConfig({
  plugins: [react()],
  test: {
    environment: "happy-dom",
    globals: true,
    setupFiles: ["./src/__tests__/setup.js"],
    coverage: {
      provider: "v8",
      reporter: ["text", "json", "html"],
    },
  },
})

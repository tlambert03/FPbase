import react from "@vitejs/plugin-react"
import { defineConfig } from "vite"

const external = ["react", "react-dom"]

export default defineConfig({
  build: {
    lib: {
      entry: "src/index.tsx",
      formats: ["es"],
    },
    rollupOptions: {
      external,
      output: {
        // Manual chunking for optimal code splitting
        manualChunks: (id) => {
          // Vendor chunks
          if (id.includes("node_modules")) {
            // Separate MUI into its own chunk
            if (id.includes("@mui") || id.includes("@emotion")) {
              return "mui-vendor"
            }
            // Highcharts in its own chunk
            if (id.includes("highcharts")) {
              return "highcharts-vendor"
            }
            // TanStack Query in its own chunk
            if (id.includes("@tanstack/react-query")) {
              return "query-vendor"
            }
            // Other vendors
            return "vendor"
          }
        },
      },
    },
    // Optimize chunk size
    chunkSizeWarningLimit: 600,
    // Enable source maps for debugging
    sourcemap: true,
    // Enable minification
    minify: "esbuild",
    target: ["es2015", "safari13"],
  },
  server: {
    proxy: {
      "/graphql": "http://127.0.0.1:8000",
      "/api": "http://127.0.0.1:8000",
    },
  },
  plugins: [react()],
  // Optimize dependencies
  optimizeDeps: {
    include: ["react", "react-dom"],
    exclude: [],
  },
})

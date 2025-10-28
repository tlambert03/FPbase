import { defineConfig } from "vite"
import react from "@vitejs/plugin-react"

const external = ["react"]

export default defineConfig({
  build: {
    lib: {
      entry: "src/index.jsx",
      formats: ["es"],
    },
    rollupOptions: { external },
  },
  server: {
    proxy: {
      "/api": "http://127.0.0.1:8000",
    },
  },
  plugins: [react()],
})

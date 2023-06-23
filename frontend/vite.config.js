import { defineConfig } from "vite"
import react from "@vitejs/plugin-react"

export default defineConfig({
  root: "./src",
  base: "/static/",
  server: {
    host: "localhost",
    port: 3000,
    open: false,
    watch: {
      usePolling: true,
      disableGlobbing: false,
    },
  },
  resolve: {
    extensions: [".js", ".jsx", ".json"],
  },
  optimizeDeps: {
    exclude: ["regenerator-runtime/runtime", "bootstrap"],
  },
  define: {
    "process.env.RESET_APP_DATA_TIMER": false
  },
  build: {
    // Relative to the root
    outDir: "../dist",
    assetsDir: "",
    manifest: true,
    emptyOutDir: true,
    target: "es2015",

    // NOTE: perhaps add images to manifest in build script... or add plugin
  // on_build_end
    plugins: [
      // with this plugin I get this problem:
      // https://stackoverflow.com/questions/75883720/504-outdated-optimize-dep-while-using-react-vite
      react(),
    ],

    rollupOptions: {
      input: {
        main: "./src/index.js",
        embedscope: "./src/embedscope.js",
        litemol: "./src/my-litemol.js",
        spectraViewer: "./src/spectra-viewer.js",
        simpleSpectraViewer: "./src/simple-spectra-viewer.js",
        microscopeForm: "./src/microscope-form.js",
        blast: "./src/blast-app.js",
      },
      output: {
        chunkFileNames: undefined,
      },
    },
  },
})

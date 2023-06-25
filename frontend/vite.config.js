import { defineConfig } from "vite"

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
    proxy: {
      "/graphql": "http://127.0.0.1:8000",
      "/api": "http://127.0.0.1:8000",
    },
  },
  resolve: {
    extensions: [".js", ".jsx", ".json"],
  },
  optimizeDeps: {
    exclude: ["regenerator-runtime/runtime", "bootstrap"],
  },
  define: {
    "process.env.RESET_APP_DATA_TIMER": false,
  },
  build: {
    // Relative to the root
    outDir: "../dist",
    assetsDir: "",
    manifest: true,
    emptyOutDir: true,
    target: "es2015",
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

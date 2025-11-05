import path from "node:path"
import { fileURLToPath } from "node:url"
import { sentryVitePlugin } from "@sentry/vite-plugin"
import react from "@vitejs/plugin-react"
import { visualizer } from "rollup-plugin-visualizer"
import { defineConfig } from "vite"
import { viteStaticCopy } from "vite-plugin-static-copy"

const __dirname = path.dirname(fileURLToPath(import.meta.url))

export default defineConfig(({ mode }) => {
  const isDev = mode === "development"
  const isTestBuild = process.env.TEST_BUILD === "1"

  return {
    // Match Django STATIC_URL
    base: "/static/",

    // Build configuration
    build: {
      // Output to frontend/dist/
      outDir: "dist",
      emptyOutDir: true,

      target: ["es2015", "safari13"],

      // Generate manifest.json for django-vite
      manifest: "manifest.json",

      // Source maps: inline for tests (easier debugging), external for production (Sentry upload)
      sourcemap: isTestBuild ? "inline" : true,

      // Minification: disabled for test builds to preserve function names and line numbers
      minify: isTestBuild ? false : "esbuild",

      // Chunk size warnings (increased for large vendor chunks)
      chunkSizeWarningLimit: 1000, // 1MB (we have large MUI/Highcharts chunks)

      // Multi-page app configuration
      rollupOptions: {
        input: {
          main: path.resolve(__dirname, "src/index.js"),
          d3Charts: path.resolve(__dirname, "src/d3-charts.js"), // D3 + chart components (lazy-loaded)
          embedscope: path.resolve(__dirname, "src/embedscope.js"),
          litemol: path.resolve(__dirname, "src/my-litemol.js"),
          spectraViewer: path.resolve(__dirname, "src/spectra-viewer.js"),
          simpleSpectraViewer: path.resolve(__dirname, "src/simple-spectra-viewer.js"),
          microscopeForm: path.resolve(__dirname, "src/microscope-form.js"),
          blast: path.resolve(__dirname, "src/blast-app.js"),
          proteinTable: path.resolve(__dirname, "src/protein-table.js"),
          scopeReport: path.resolve(__dirname, "src/scope-report.js"),
          fret: path.resolve(__dirname, "src/fret.js"),
        },

        // Let Vite handle code splitting automatically
        // Manual chunking was causing initialization order issues
        output: {},
      },
    },

    // Development server
    server: {
      port: 5173,
      strictPort: true,
      origin: "http://localhost:5173",
      cors: true,
      hmr: {
        protocol: "ws",
        host: "localhost",
      },
    },

    resolve: {
      alias: {
        "@fpbase/spectra": path.resolve(__dirname, "../packages/spectra/src/index.tsx"),
        "@fpbase/blast": path.resolve(__dirname, "../packages/blast/src/index.js"),
        "@fpbase/protein-table": path.resolve(__dirname, "../packages/protein-table/src/index.jsx"),
        // jQuery loaded from CDN - no alias needed
      },
    },

    // Plugins
    plugins: [
      // React with Fast Refresh
      react({
        // Include .js files for JSX processing (not just .jsx)
        include: /\.(jsx|js|tsx|ts)$/,
      }),

      // Copy static files (microscope.js)
      viteStaticCopy({
        targets: [
          {
            src: "../backend/fpbase/static/js/microscope.js",
            dest: "js",
          },
        ],
      }),

      // Sentry source map upload (production only, and only if auth token is set)
      !isDev &&
        process.env.SENTRY_AUTH_TOKEN &&
        sentryVitePlugin({
          org: "talley-lambert",
          project: "fpbase",
          authToken: process.env.SENTRY_AUTH_TOKEN,
          release: process.env.HEROKU_SLUG_COMMIT,
          telemetry: false,

          // Tree-shake optional Sentry code to reduce bundle size
          bundleSizeOptimizations: {
            excludeDebugStatements: true, // Remove debug logging (dev-only code)
            excludePerformanceMonitoring: false, // Keep - we use browserTracingIntegration
            excludeReplayIframe: true, // Remove - we don't record iframes
            excludeReplayShadowDom: true, // Remove - we don't use shadow DOM
            excludeReplayWorker: true, // Use main thread compression (smaller bundle)
          },
        }),

      // Bundle analyzer (production only, use `ANALYZE=1 pnpm build` to generate stats.html)
      !isDev &&
        process.env.ANALYZE &&
        visualizer({
          filename: "stats.html",
          open: true,
          gzipSize: true,
          brotliSize: true,
        }),
    ].filter(Boolean),

    // CSS configuration
    css: {
      postcss: {
        plugins: [require("autoprefixer"), require("cssnano")],
      },
      preprocessorOptions: {
        scss: {
          silenceDeprecations: ["import", "global-builtin", "color-functions", "abs-percent"],
        },
      },
    },

    // Define environment variables
    define: {
      "process.env.NODE_ENV": JSON.stringify(mode),
      "process.env.SENTRY_DSN": JSON.stringify(process.env.SENTRY_DSN || ""),
      "process.env.HEROKU_SLUG_COMMIT": JSON.stringify(process.env.HEROKU_SLUG_COMMIT || ""),
    },

    // Configure esbuild to handle JSX/TSX in .js, .jsx, .ts, and .tsx files
    esbuild: {
      loader: "tsx",
      include: /src\/.*\.[jt]sx?$/,
      exclude: [],
    },

    // Optimize dependencies
    optimizeDeps: {
      // Exclude jQuery - loaded from CDN
      exclude: ["jquery"],
      include: ["process/browser"],
      esbuildOptions: {
        // Handle JSX/TSX in .js and .ts files during dependency scanning
        loader: {
          ".js": "jsx",
          ".ts": "tsx",
        },
        define: {
          global: "globalThis",
        },
      },
    },
  }
})

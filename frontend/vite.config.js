import path from "node:path"
import { fileURLToPath } from "node:url"
import { sentryVitePlugin } from "@sentry/vite-plugin"
import react from "@vitejs/plugin-react"
import { defineConfig } from "vite"
import { viteStaticCopy } from "vite-plugin-static-copy"

const __dirname = path.dirname(fileURLToPath(import.meta.url))

export default defineConfig(({ mode }) => {
  const isDev = mode === "development"

  return {
    // Match Django STATIC_URL
    base: "/static/",

    // Build configuration
    build: {
      // Output to frontend/dist/
      outDir: "dist",
      emptyOutDir: true,

      // Generate manifest.json for django-vite
      manifest: "manifest.json",

      // Source maps for Sentry
      sourcemap: true,

      // Multi-page app configuration
      rollupOptions: {
        input: {
          main: path.resolve(__dirname, "src/index.js"),
          embedscope: path.resolve(__dirname, "src/embedscope.js"),
          litemol: path.resolve(__dirname, "src/my-litemol.js"),
          spectraViewer: path.resolve(__dirname, "src/spectra-viewer.js"),
          simpleSpectraViewer: path.resolve(__dirname, "src/simple-spectra-viewer.js"),
          microscopeForm: path.resolve(__dirname, "src/microscope-form.js"),
          blast: path.resolve(__dirname, "src/blast-app.js"),
          proteinTable: path.resolve(__dirname, "src/protein-table.js"),
        },

        // Manual code splitting (hybrid approach)
        output: {
          manualChunks(id) {
            // Sentry (shared across all)
            if (id.includes("node_modules/@sentry")) {
              return "vendor-sentry"
            }

            // React (shared, but NOT in embedscope)
            if (id.includes("node_modules/react") || id.includes("node_modules/scheduler")) {
              return "vendor-react"
            }

            // jQuery (shared across multiple)
            if (id.includes("node_modules/jquery")) {
              return "vendor-jquery"
            }

            // D3 v7 (EXCLUDE from embedscope - it uses CDN D3 v3)
            if (id.includes("node_modules/d3")) {
              return "vendor-d3"
            }
          },
        },
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

    // Resolve aliases (match webpack)
    resolve: {
      alias: {
        "@fpbase/spectra": path.resolve(__dirname, "../packages/spectra/src/index.jsx"),
        "@fpbase/blast": path.resolve(__dirname, "../packages/blast/src/index.js"),
        "@fpbase/protein-table": path.resolve(__dirname, "../packages/protein-table/src/index.jsx"),
        jquery: path.resolve(__dirname, "node_modules/jquery/src/jquery"),
      },
    },

    // Plugins
    plugins: [
      // React with Fast Refresh
      react({
        // Include .js files for JSX processing (not just .jsx)
        include: /\.(jsx|js|tsx|ts)$/,
        babel: {
          plugins: ["@babel/plugin-syntax-dynamic-import"],
        },
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

      // Sentry source map upload (production only)
      !isDev &&
        sentryVitePlugin({
          org: "talley-lambert",
          project: "fpbase",
          authToken: process.env.SENTRY_AUTH_TOKEN,
          release: process.env.HEROKU_SLUG_COMMIT,
        }),
    ].filter(Boolean),

    // CSS configuration
    css: {
      postcss: {
        plugins: [require("autoprefixer"), require("cssnano")],
      },
      preprocessorOptions: {
        scss: {
          // Silence deprecation warnings similar to webpack config
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

    // Configure esbuild to handle JSX in .js and .jsx files
    esbuild: {
      loader: "jsx",
      include: /src\/.*\.[jt]sx?$/,
      exclude: [],
    },

    // Optimize dependencies
    optimizeDeps: {
      include: ["jquery", "process/browser"],
      esbuildOptions: {
        // Inject global jQuery (match webpack ProvidePlugin)
        define: {
          global: "globalThis",
        },
      },
    },
  }
})

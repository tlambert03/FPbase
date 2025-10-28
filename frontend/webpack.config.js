const path = require("path")
const webpack = require("webpack")
const BundleTracker = require("webpack-bundle-tracker")
const { CleanWebpackPlugin } = require("clean-webpack-plugin")
const CopyPlugin = require("copy-webpack-plugin")
const MiniCssExtractPlugin = require("mini-css-extract-plugin")
const CssMinimizerPlugin = require("css-minimizer-webpack-plugin")
const TerserJSPlugin = require("terser-webpack-plugin")
const { sentryWebpackPlugin } = require("@sentry/webpack-plugin")

const devMode = process.env.NODE_ENV !== "production"
const hotReload = process.env.HOT_RELOAD === "1"

const styleRule = {
  test: /\.(sa|sc|c)ss$/,
  use: [
    MiniCssExtractPlugin.loader,
    {
      loader: "css-loader",
      options: { sourceMap: true },
    },
    {
      loader: "postcss-loader",
      options: {
        postcssOptions: {
          plugins: ["autoprefixer", "cssnano"],
        },
      },
    },
    {
      loader: "sass-loader",
      options: {
        sassOptions: {
          silenceDeprecations: ['import', 'global-builtin', 'color-functions', 'abs-percent'],
        },
      },
    },
  ],
}

const jsRule = {
  test: /\.jsx?$/,
  exclude: /node_modules\/(?!(qs|debounce-raf)\/).*/,
  use: {
    loader: "babel-loader",
    options: {
      presets: [
        [
          "@babel/preset-env",
          {
            useBuiltIns: "entry",
            corejs: 3,
          },
        ],
        [
          "@babel/preset-react",
          {
            runtime: "automatic",
          },
        ],
      ],
      plugins: ["@babel/plugin-syntax-dynamic-import"],
    },
  },
}

const assetRule = {
  test: /.(jpe?g|png|woff(2)?|eot|ttf|svg)$/,
  type: "asset",
}

const plugins = [
  new webpack.ProgressPlugin(),
  new webpack.ProvidePlugin({
    $: "jquery",
    jQuery: "jquery",
    process: "process/browser",
  }),
  // Inject environment variables into the bundle for runtime access
  new webpack.DefinePlugin({
    'process.env.NODE_ENV': JSON.stringify(process.env.NODE_ENV || 'development'),
    'process.env.SENTRY_DSN': JSON.stringify(process.env.SENTRY_DSN || ''),
    'process.env.HEROKU_SLUG_COMMIT': JSON.stringify(process.env.HEROKU_SLUG_COMMIT || ''),
  }),
  // new webpack.IgnorePlugin(/vertx/),
  new BundleTracker({
    filename: "webpack-stats.json",
  }),
  new MiniCssExtractPlugin({
    filename: "[name].[contenthash].css",
    chunkFilename: "[id].[chunkhash].css",
  }),
  new CleanWebpackPlugin(),
  new CopyPlugin({
    patterns: [
      {
        from: path.resolve(__dirname, "../backend/fpbase/static/js/microscope.js"),
        to: "js/microscope.js",
      },
    ],
  }),
]

if (!devMode) {
  plugins.push(
    sentryWebpackPlugin({
      org: "talley-lambert",
      project: "fpbase",
      authToken: process.env.SENTRY_AUTH_TOKEN,
      include: "./dist",
      ignore: ["node_modules", "webpack.config.js"],
      release: process.env.HEROKU_SLUG_COMMIT,
    })
  )
}

module.exports = {
  context: __dirname,
  entry: {
    main: "./src/index.js",
    embedscope: "./src/embedscope.js",
    litemol: "./src/my-litemol.js",
    spectraViewer: "./src/spectra-viewer.js",
    simpleSpectraViewer: "./src/simple-spectra-viewer.js",
    microscopeForm: "./src/microscope-form.js",
    blast: "./src/blast-app.js",
    proteinTable: "./src/protein-table.js",
  },
  output: {
    path: path.resolve("./dist/"),
    filename: devMode ? "[name].js" : "[name].[contenthash].js",
    publicPath: hotReload ? "http://localhost:8080/static/" : "/static/",
    chunkFilename: devMode ? "[name].js" : "[name].[chunkhash].js",
  },
  resolve: {
    extensions: [".webpack.js", ".web.js", ".mjs", ".js", ".jsx", ".json"],
    alias: {
      jquery: "jquery/src/jquery",
      "@fpbase/spectra": path.resolve(
        __dirname,
        "../packages/spectra/src/index.jsx"
      ),
      "@fpbase/blast": path.resolve(
        __dirname,
        "../packages/blast/src/index.js"
      ),
      "@fpbase/protein-table": path.resolve(
        __dirname,
        "../packages/protein-table/src/index.jsx"
      ),
      "@emotion/react": require.resolve("@emotion/react"),
      "@emotion/styled": require.resolve("@emotion/styled"),
    },
    fallback: {
      url: require.resolve("url/"),
    },
  },
  devtool: devMode ? "eval-cheap-module-source-map" : "source-map",
  devServer: {
    port: 8080,
    headers: {
      "Access-Control-Allow-Origin": "*",
    },
    static: {
      directory: path.resolve(__dirname, "../backend/fpbase/static"),
      publicPath: "/",
    },
  },
  module: {
    rules: [
      {
        test: /\.m?js$/,
        resolve: { fullySpecified: false },
        include: /node_modules/,
      },
      jsRule,
      styleRule,
      assetRule,
    ],
  },
  plugins,
  optimization: {
    // Share webpack runtime across all bundles
    runtimeChunk: "single",

    minimizer: [
      new TerserJSPlugin(),
      new CssMinimizerPlugin({
        minimizerOptions: {
          preset: ["default", { svgo: false }],
        },
      }),
    ],

    // Improved code splitting for shared dependencies
    splitChunks: {
      chunks: "all",
      cacheGroups: {
        // Extract Sentry into its own chunk (shared across all bundles)
        sentry: {
          test: /[\\/]node_modules[\\/]@sentry[\\/]/,
          name: "sentry",
          priority: 30,
          reuseExistingChunk: true,
        },
        // Extract React and React-DOM (used by multiple bundles)
        react: {
          test: /[\\/]node_modules[\\/](react|react-dom|scheduler)[\\/]/,
          name: "react-vendor",
          priority: 25,
          reuseExistingChunk: true,
        },
        // Extract jQuery (used by multiple bundles)
        jquery: {
          test: /[\\/]node_modules[\\/]jquery[\\/]/,
          name: "jquery",
          priority: 20,
          reuseExistingChunk: true,
        },
        // Extract D3 separately - NOT used by embedscope (uses CDN D3 v3 instead)
        d3: {
          test: /[\\/]node_modules[\\/]d3.*[\\/]/,
          name: "d3-vendor",
          priority: 15,
          chunks: (chunk) => chunk.name !== "embedscope",
          reuseExistingChunk: true,
        },
        // Extract other large vendor libraries (excluding D3 which has its own chunk)
        vendors: {
          test(module) {
            // Exclude D3 from vendors chunk - it has its own d3-vendor chunk
            return (
              module.resource &&
              /[\\/]node_modules[\\/]/.test(module.resource) &&
              !/[\\/]node_modules[\\/]d3.*[\\/]/.test(module.resource)
            )
          },
          name: "vendors",
          priority: 10,
          reuseExistingChunk: true,
        },
        // Extract common code shared between entry points
        common: {
          minChunks: 2,
          name: "common",
          priority: 5,
          reuseExistingChunk: true,
        },
      },
    },
  },
}

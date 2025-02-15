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
    "sass-loader",
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
        "@babel/preset-react",
      ],
      plugins: ["@babel/plugin-syntax-dynamic-import"],
    },
  },
}

const assetRule = {
  test: /.(jpe?g|png|woff(2)?|eot|ttf|svg)$/,
  loader: "file-loader",
  type: "asset",
}

const plugins = [
  new webpack.ProgressPlugin(),
  new webpack.ProvidePlugin({
    $: "jquery",
    jQuery: "jquery",
    process: "process/browser",
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
]

if (!devMode) {
  plugins.push(
    sentryWebpackPlugin({
      org: "talley-lambert",
      project: "fpbase",
      authToken: process.env.SENTRY_AUTH_TOKEN,
      include: "./dist",
      ignore: ["node_modules", "webpack.config.js"],
      release: process.env.SOURCE_VERSION,
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
    minimizer: [new TerserJSPlugin(), new CssMinimizerPlugin()],
    splitChunks: {
      cacheGroups: {
        commons: {
          test: /[\\/]node_modules[\\/]/,
          name: "vendor",
          chunks: "initial",
        },
      },
    },
  },
}

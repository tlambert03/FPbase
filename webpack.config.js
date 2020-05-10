const path = require("path")
const webpack = require("webpack")
const autoprefixer = require("autoprefixer")
const { BundleAnalyzerPlugin } = require("webpack-bundle-analyzer")
const BundleTracker = require("webpack-bundle-tracker")
const { CleanWebpackPlugin } = require("clean-webpack-plugin")
const CopyWebpackPlugin = require("copy-webpack-plugin")
const ImageminWebpackPlugin = require("imagemin-webpack-plugin").default
const ImageminWebP = require("imagemin-webp")
const MiniCssExtractPlugin = require("mini-css-extract-plugin")
const OptimizeCSSAssetsPlugin = require("optimize-css-assets-webpack-plugin")
const SentryCliPlugin = require("@sentry/webpack-plugin")
const TerserJSPlugin = require("terser-webpack-plugin")
const CSSnano = require("cssnano")

const devMode = process.env.NODE_ENV !== "production"
const hotReload = process.env.HOT_RELOAD === "1"

const styleRule = {
  test: /\.(sa|sc|c)ss$/,
  use: [
    MiniCssExtractPlugin.loader,
    {
      loader: "css-loader",
      options: {
        sourceMap: true,
      },
    },
    {
      loader: "postcss-loader",
      options: {
        plugins: () => [autoprefixer(), CSSnano],
      },
    },
    "sass-loader",
  ],
}

const jsRule = {
  test: /\.jsx?$/,
  exclude: /node_modules\/(?!(smiles-drawer|qs|debounce-raf)\/).*/,
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
  test: /.(jpg|png|woff(2)?|eot|ttf|svg)$/,
  loader: "file-loader",
}

const plugins = [
  new webpack.ProgressPlugin(),
  new webpack.ProvidePlugin({
    $: "jquery",
    jQuery: "jquery",
  }),
  new webpack.IgnorePlugin(/vertx/),
  new BundleTracker({
    filename: "./webpack-stats.json",
  }),
  new MiniCssExtractPlugin({
    filename: "[name].[contenthash].css",
    chunkFilename: "[id].[chunkhash].css",
  }),
  new BundleAnalyzerPlugin({
    analyzerMode: "static",
    openAnalyzer: false,
  }),
  new CleanWebpackPlugin(),
  new CopyWebpackPlugin([
    {
      from: "./static/src/images/**/*",
      to: path.resolve("./static/dist/images/[name].webp"),
      toType: "template",
    },
  ]),
  new CopyWebpackPlugin([
    {
      from: "./static/src/js/sentry.*.js",
      to: path.resolve("./static/dist/sentry.js"),
    },
  ]),
  new ImageminWebpackPlugin({
    test: /\.(webp)$/i,
    plugins: [
      ImageminWebP({
        quality: 90,
        sharpness: 1,
      }),
    ],
  }),
  new CopyWebpackPlugin([
    {
      from: "./static/src/images/**/*",
      to: path.resolve("./static/dist/images/[name].[ext]"),
      toType: "template",
    },
  ]),
]

if (devMode) {
  plugins.push(new webpack.HotModuleReplacementPlugin())
}

if (devMode) {
  styleRule.use = ["css-hot-loader", ...styleRule.use]
} else {
  plugins.push(
    new webpack.EnvironmentPlugin({
      NODE_ENV: "development",
      SOURCE_VERSION: false,
      SENTRY_DSN: false,
      SENTRY_AUTH_TOKEN: false,
    })
  )

  if (process.env.SENTRY_AUTH_TOKEN && !process.env.CI) {
    plugins.push(
      new SentryCliPlugin({
        include: "static/",
        release: process.env.SOURCE_VERSION,
        ignore: [
          "node_modules",
          "webpack.config.js",
          "static/src/js/pdb/LiteMol-plugin.js",
        ],
      })
    )
  }
}

module.exports = {
  context: __dirname,
  entry: {
    main: "./static/src/index.js",
    embedscope: "./static/src/embedscope.js",
    litemol: "./static/src/my-litemol.js",
    spectraViewer: "./static/src/spectra-viewer.js",
    simpleSpectraViewer: "./static/src/simple-spectra-viewer.js",
    microscopeForm: "./static/src/microscope-form.js",
    blast: "./static/src/blast-app.js",
  },
  output: {
    path: path.resolve("./static/dist/"),
    filename: devMode ? "[name].js" : "[name].[hash].js",
    publicPath: hotReload ? "http://localhost:8080/static/" : "/static/",
    // publicPath: hotReload ? 'http://10.0.2.2:8080/static/' : '/static/',
    chunkFilename: devMode ? "[name].js" : "[name].[chunkhash].js",
  },
  resolve: {
    extensions: [".js", ".jsx"],
    alias: {
      jquery: "jquery/src/jquery",
    },
  },
  devtool: devMode ? "cheap-module-eval-source-map" : "source-map",
  devServer: {
    hot: true,
    overlay: true,
    quiet: false,
    port: 8080,
    headers: {
      "Access-Control-Allow-Origin": "*",
    },
  },
  module: {
    rules: [jsRule, styleRule, assetRule],
  },
  externals: {
    Sentry: "Sentry",
  },
  plugins,
  optimization: {
    // minimizer: [
    //   new UglifyJsPlugin({
    //     cache: true,
    //     parallel: true,
    //     sourceMap: true // set to true if you want JS source maps
    //   }),
    //   new OptimizeCSSAssetsPlugin({})
    // ],
    minimizer: [
      new TerserJSPlugin({ cache: true, parallel: true, sourceMap: true }),
      new OptimizeCSSAssetsPlugin({}),
    ],
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

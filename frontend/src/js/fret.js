import Highcharts from "highcharts"
import "highcharts/modules/no-data-to-display"
import { fetchWithSentry } from "./ajax-sentry"
import { icon } from "./icons.js"

const $ = window.jQuery // jQuery loaded from CDN

export default function initFRET() {
  var chart
  var data = []
  var localData = {}
  var donorEM
  var acceptorEX
  var options = {
    showArea: true,
    minwave: 350,
    maxwave: 750,
    startingBrush: [350, 750],
    autoscaleBrush: false,
    exNormWave: undefined,
    scale: "linear",
    hide2p: true,
    scaleToEC: false,
    scaleToQY: false,
  }

  //initialize chart with Highcharts
  chart = Highcharts.chart("spectra", {
    chart: {
      type: "line",
      animation: { duration: 300 },
      backgroundColor: "#ffffff",
      height: 350,
    },
    title: { text: null },
    xAxis: {
      title: { text: "Wavelength (nm)" },
      min: 350,
      max: 750,
      tickLength: 0,
      gridLineWidth: 0,
    },
    yAxis: {
      title: { text: null }, // Hide Y axis label
      min: 0,
      max: 1,
      labels: {
        formatter: function () {
          return `${Math.round(this.value * 100)}%`
        },
      },
      gridLineWidth: 0, // Remove horizontal grid lines
    },
    tooltip: {
      shared: !window.mobilecheck(),
      valueDecimals: 1,
      valueSuffix: "%",
      valueFormatter: (value) => {
        if (value !== null && value !== undefined) {
          return `${Math.round(value * 1000) / 10}%`
        } else {
          return "--"
        }
      },
    },
    legend: {
      enabled: true,
      align: "right",
      verticalAlign: "top",
      layout: "horizontal",
      x: 0,
      y: 0,
      floating: false,
      symbolWidth: 12,
      symbolHeight: 12,
      symbolRadius: 12,
    },
    plotOptions: {
      series: {
        legendSymbol: "rectangle",
      },
      line: {
        animation: false,
        marker: { enabled: false },
        lineWidth: 2,
      },
      area: {
        animation: false,
        marker: { enabled: false },
        lineWidth: 2,
      },
    },
    series: [],
    lang: {
      noData: "Choose a donor and acceptor below",
    },
    noData: {
      style: {
        fontWeight: "bold",
        fontSize: "15px",
        color: "#303030",
      },
    },
    credits: { enabled: false },
    accessibility: { enabled: false },
  })

  // Add resize handler
  $(window).on("resize", () => {
    if (chart) {
      chart.reflow()
    }
  })

  $("#donor-select").select2({
    theme: "bootstrap-5",
    containerCssClass: ":all:",
    width: "resolve",
  })
  $("#acceptor-select").select2({
    theme: "bootstrap-5",
    containerCssClass: ":all:",
    width: "resolve",
  })

  function getData(slug) {
    var dfd = $.Deferred()
    // download if not already downloaded
    if (!(slug in localData)) {
      fetchWithSentry(`/spectra/${slug}`, {
        // Legacy header required by Django is_ajax() check in dual-purpose endpoints
        headers: { "X-Requested-With": "XMLHttpRequest" },
      })
        .then((response) => {
          if (!response.ok) {
            throw new Error(`HTTP ${response.status}: ${response.statusText}`)
          }
          return response.json()
        })
        .then((d) => {
          for (let n = 0; n < d.spectra.length; n++) {
            if (d.spectra[n].type !== "2p") {
              d.spectra[n] = padDataLimits(d.spectra[n])
            }
          }
          localData[slug] = d.spectra
          dfd.resolve(localData[slug])
        })
        .catch((error) => {
          console.error(`Failed to get spectra for ${slug}:`, error)
          dfd.reject(0)
        })
    } else {
      // otherwise pull from local dict
      dfd.resolve(localData[slug])
    }
    return dfd.promise()
  }

  function padDataLimits(d) {
    for (let i = d.minwave - 1; i >= options.minwave; i--) {
      d.values.unshift({ x: i, y: 0 })
    }
    for (let n = d.maxwave + 1; n <= Math.max(options.maxwave, 1000); n++) {
      d.values.push({ x: n, y: 0 })
    }
    return d
  }

  function dataHasKey(key) {
    return $.grep(data, (obj) => obj.key === key).length > 0
  }

  function dataItemMatching(filter, d) {
    d = d || data
    return d.filter((item) => {
      for (var key in filter) {
        if (typeof item[key] === "undefined" || item[key] !== filter[key]) return false
      }
      return true
    })
  }

  function addItem(slug, subtype) {
    if (slug === "") {
      return $.Deferred((def) => {
        def.resolve()
      }).promise()
    }
    // add spectral data to array
    subtype = subtype || false

    return getData(slug)
      .then((d) => {
        for (let i = 0; i < d.length; i++) {
          if ((d[i].type !== "2p") & !dataHasKey(d[i].key)) {
            data.push(JSON.parse(JSON.stringify(d[i]))) // make a copy of the object
          }
        }
      })
      .fail((_d) => {
        console.log("item not found")
      })
  }

  function spectral_product(ar1, ar2) {
    // calculate product of two spectra.values
    var output = []
    //var step = ar1[1].x - ar1[0].x;
    var left = Math.max(ar1[0].x, ar2[0].x)
    var right = Math.min(ar1[ar1.length - 1].x, ar2[ar2.length - 1].x)

    var a1 = ar1.slice(
      ar1.findIndex((i) => i.x === left),
      ar1.findIndex((i) => i.x === right)
    )
    var a2 = ar2.slice(
      ar2.findIndex((i) => i.x === left),
      ar2.findIndex((i) => i.x === right)
    )
    for (let i = 0; i < a1.length; i++) {
      output.push({ x: a1[i].x, y: a1[i].y * a2[i].y })
    }

    return output
  }

  function forster_distance(ar1, ar2, k2, ni) {
    k2 = k2 || 2 / 3
    ni = ni || 1.33

    var donorQY = ar1.scalar
    var accECmax = ar2.scalar
    var a1wavemap = ar1.values.reduce((acc, cur) => {
      acc[cur.x] = cur.y
      return acc
    }, {})
    var a2wavemap = ar2.values.reduce((acc, cur) => {
      acc[cur.x] = cur.y
      return acc
    }, {})
    var donorsum = ar1.values.reduce((a, b) => a + b.y, 0)
    var startingwave = Math.max(ar1.minwave, ar2.minwave)
    var endingwave = Math.min(ar1.maxwave, ar2.maxwave)
    var step = ar1.values[1].x - ar1.values[0].x
    let overlapIntgrl = 0
    for (let wave = startingwave; wave <= endingwave; wave += step) {
      overlapIntgrl +=
        ((wave * 1e-7) ** 4 * a1wavemap[wave] * accECmax * a2wavemap[wave]) / donorsum
    }
    var r = (8.8 * 1e-25 * k2 * donorQY * ni ** -4 * overlapIntgrl) ** (1 / 6) * 1e7
    return [overlapIntgrl, r]
  }

  $("#ni-input").on("change", function () {
    var input = $(this)
    var val = input.val()
    input.val(Math.max(Math.min(val, 1.6), 0.5))
    updateTable()
  })

  $("#k2-input").on("change", function () {
    var input = $(this)
    var val = input.val()
    input.val(Math.max(Math.min(val, 4), 0))
    updateTable()
  })

  function updateTable() {
    var r = forster_distance(donorEM, acceptorEX, $("#k2-input").val(), $("#ni-input").val())
    $("#overlapIntgrl").text(Math.round(r[0] * 1e15) / 100)
    $("#r0").text(Math.round(r[1] * 1000) / 100)
    $("#r0QYA").text(Math.round(r[1] * $("#QYA").text() * 1000) / 100)
  }

  function updateChart() {
    // Convert data to Highcharts format
    var series = data.map((item) => {
      var isFadedFret = item.classed?.includes("faded-fret")
      var isFretOverlap = item.classed?.includes("fret-overlap")

      var seriesConfig = {
        name: item.key,
        data: item.values.map((v) => [v.x, v.y]),
        type: item.area ? "area" : "line",
        className: item.classed,
        zIndex: isFretOverlap ? 10 : isFadedFret ? 1 : 5,
      }

      // Apply faded-fret styling
      if (isFadedFret) {
        seriesConfig.dashStyle = "Dash" // Highcharts dash style (5,5)
        seriesConfig.lineColor = "rgba(0, 0, 0, 0.2)" // Gray dashed line
        seriesConfig.lineWidth = 1.2
        if (item.area) {
          seriesConfig.fillOpacity = 0.15 // Keep original color but very transparent
        }
        // Keep the original item.color for the fill
        if (item.color && typeof item.color === "string" && !item.color.startsWith("url(")) {
          seriesConfig.color = item.color
        }
      } else {
        // Non-faded items: no line, just filled area
        seriesConfig.lineWidth = 0
      }

      // Handle color - could be a string, object (pattern), or undefined
      if (item.color && !isFadedFret) {
        if (typeof item.color === "string" && item.color.startsWith("url(")) {
          // Skip old-style SVG pattern references
          seriesConfig.color = null
        } else {
          // Use the color directly (string or pattern object)
          seriesConfig.color = item.color
        }
      }

      // Set fillOpacity for area charts (but not for pattern fills or faded-fret)
      if (item.area && !isFadedFret && (!item.color || typeof item.color === "string")) {
        seriesConfig.fillOpacity = 0.5
      }

      return seriesConfig
    })

    // Update chart with new series
    while (chart.series.length > 0) {
      chart.series[0].remove(false)
    }
    series.forEach((s) => {
      chart.addSeries(s, false)
    })
    chart.redraw()
  }

  // main function when data-selector has been changed
  $(".data-selector").change((_event) => {
    data.splice(0, 10)
    var donorslug = $("#donor-select :selected").val()
    var acceptorslug = $("#acceptor-select :selected").val()
    $.when(addItem(donorslug), addItem(acceptorslug)).then(() => {
      if (donorslug || acceptorslug) {
        $(".table-wrapper").removeAttr("hidden").show()
      }
      if (donorslug) {
        donorEM = dataItemMatching({ slug: donorslug, type: "em" })[0]
        donorEM.classed = `${donorEM.classed || ""} faded-fret`
        $("#QYD").text(donorEM.scalar)
      } else {
        $("#QYD").text("")
      }
      if (acceptorslug) {
        acceptorEX = dataItemMatching({ slug: acceptorslug, type: "ex" })[0]
        acceptorEX.classed = `${acceptorEX.classed || ""} faded-fret`
        $("#ECA").text(acceptorEX.scalar.toLocaleString())

        const acceptorEM = dataItemMatching({ slug: acceptorslug, type: "em" })[0]
        $("#QYA").text(acceptorEM ? acceptorEM.scalar : "")
      } else {
        $("#ECA, #QYA").text("")
      }
      if (donorslug && acceptorslug) {
        // Dynamically import pattern-fill module only when needed
        import("highcharts/modules/pattern-fill").then(() => {
          data.push({
            key: "Overlap",
            values: spectral_product(donorEM.values, acceptorEX.values),
            area: true,
            color: {
              pattern: {
                path: {
                  d: "M -1 1 l 2 -2 M 0 10 l 10 -10 M 9 11 l 2 -2",
                  stroke: "white",
                  strokeWidth: 3,
                },
                width: 10,
                height: 10,
                backgroundColor: "#00000092",
                opacity: 1,
              },
            },
            classed: "fret-overlap",
            type: "overlap",
          })
          updateChart()
        })
        updateTable(donorEM, acceptorEX)
      } else {
        $("#overlapIntgrl, #r0").text("")
        updateChart()
      }
    })
  })

  $("body").on("click", ".load-button", function () {
    var donorslug = `${$(this).closest("tr").find("td:nth-child(2) a").attr("href").split("/")[2]}_default`
    var accslug = `${$(this).closest("tr").find("td:nth-child(3) a").attr("href").split("/")[2]}_default`
    $("#donor-select").val(donorslug).trigger("change.select2")
    $("#acceptor-select").val(accslug).change()
  })

  /* Custom filtering function which will search data in column four between two values */
  $.fn.dataTable.ext.search.push((_settings, data, _dataIndex) => {
    var min = parseFloat($("#minQYAinput").val())
    var minlam = parseInt($("#minLambdaSep").val(), 10)
    var qya = parseFloat(data[7]) || 0 // use data for the age column
    var lambsep = parseFloat(data[4]) || 0 // use data for the age column

    if (min > 1) {
      min = 1
      $("#minQYAinput").val(1)
    }
    if (min < 0) {
      min = 0
      $("#minQYAinput").val(0)
    }

    if ((Number.isNaN(min) || qya >= min) && (Number.isNaN(minlam) || lambsep >= minlam)) {
      return true
    }
    return false
  })

  $(document).ready(() => {
    var getData = (data, callback, settings) => {
      fetchWithSentry("", {
        // Legacy header required by Django is_ajax() check in dual-purpose endpoints
        headers: { "X-Requested-With": "XMLHttpRequest" },
      })
        .then((response) => response.json())
        .then((d) => {
          if (d.data === null) {
            setTimeout(getData.bind(this, data, callback, settings), 1500)
          } else {
            callback({
              data: d.data,
              recordsTotal: d.data.length,
              recordsFiltered: d.data.length,
            })
          }
        })
        .catch((error) => {
          console.error("Failed to get FRET data:", error)
        })
    }

    var fretTable = $("#fret_report").DataTable({
      ajax: getData,
      buttons: ["copy", "excel", "pdf"],
      scrollX: true,
      pageLength: 25,
      order: [[10, "desc"]],
      language: {
        emptyTable: "No Data received from server...",
        loadingRecords: "Recalculating FRET efficiencies across database...  Please wait.",
      },
      update: function () {
        this._positions()
        this._scroll(false)
      },
      columns: [
        {
          data: () =>
            `<button class="btn btn-sm btn-outline bg-transparent load-button">${icon("view", "text-secondary")}</button>`,
          width: "1px",
          orderable: false,
        },
        { data: "donor", width: "20px" },
        { data: "acceptor", width: "20px" },
        { data: "donorPeak", width: "10px" },
        { data: "emdist", width: "10px" },
        { data: "donorQY" },
        { data: "acceptorEC" },
        { data: "acceptorQY" },
        { data: "overlap" },
        { data: "forster" },
        { data: "forsterQYA" },
      ],
    })

    var sel = $("#fret_report_wrapper .row:first-child div.col-md-6")
    sel.removeClass("col-md-6").addClass("col-md-4")

    var e = $("#fret_report_wrapper .row:first-child div.col-md-4:first-child")
    e.removeClass("col-md-4").addClass("col-md-3")
    var D1 = $("<div>", {
      class: "col-md-2 col-sm-6 col-xs-12",
      style: "margin-top: -2px;",
    })
      .append(
        $("<div>", {
          class: "input-group input-group-sm d-flex justify-content-center pb-3",
        })
          .append(
            $("<div>", { class: "" }).append(
              $("<span>", { class: "input-group-text", id: "minQYA" }).html("min QY<sub>A</sub>")
            )
          )
          .append(
            $("<input>", {
              type: "number",
              step: 0.01,
              min: 0,
              max: 1,
              class: "form-control",
              "aria-describedby": "minQYA",
              value: 0.4,
              id: "minQYAinput",
            })
          )
      )
      .insertAfter(e)

    $("<div>", {
      class: "col-md-3 col-sm-6 col-xs-12",
      style: "margin-top: -2px;",
    })
      .append(
        $("<div>", {
          class: "input-group input-group-sm d-flex justify-content-center pb-3",
        })
          .append(
            $("<div>", { class: "" }).append(
              $("<span>", {
                class: "input-group-text",
                id: "minLambdaSepLab",
              }).html('min &Delta;&lambda;<sub class="small">em</sub>')
            )
          )
          .append(
            $("<input>", {
              type: "number",
              min: -200,
              max: 200,
              class: "form-control",
              "aria-describedby": "minLambdaSepLab",
              value: 20,
              id: "minLambdaSep",
            })
          )
      )
      .insertAfter(D1)

    $("#minQYAinput, #minLambdaSep").keyup(() => {
      fretTable.draw()
    })
  })

  var topofDiv = $(".spectra-wrapper").offset().top
  $(window).scroll(() => {
    if ($(window).scrollTop() > topofDiv) {
      $(".spectra-wrapper").css("box-shadow", "0 12px 8px -8px rgba(0,0,0,.2)")
    } else {
      $(".spectra-wrapper").css("box-shadow", "none")
    }
  })
}

import Highcharts from "highcharts"

const $ = window.jQuery // jQuery loaded from CDN

import ProgressBar from "progressbar.js"

function compare(a, b) {
  if (a.key < b.key) return 1
  if (a.key > b.key) return -1
  return 0
}

function titleCase(str) {
  str = str.replace(/_/g, " ").toLowerCase().split(" ")
  for (let i = 0; i < str.length; i++) {
    str[i] = str[i].charAt(0).toUpperCase() + str[i].slice(1)
  }
  return str.join(" ")
}

function slugify(text) {
  return text
    .toString()
    .toLowerCase()
    .replace(/\s+/g, "-") // Replace spaces with -
    .replace(/[^\w-]+/g, "") // Remove all non-word chars
    .replace(/--+/g, "-") // Replace multiple - with single -
    .replace(/^-+/, "") // Trim - from start of text
    .replace(/-+$/, "") // Trim - from end of text
}

$.fn.extend({
  scopeReport: (config) => {
    var CSRF_TOKEN = config.csrfToken
    var SCOPE_ID = config.scopeID
    var SCOPE_URL = config.scopeURL
    var JOB_ID = null
    var OUTDATED
    if (config.outdated.length) {
      OUTDATED = config.outdated
    }
    var interval = config.interval || 3500
    var FLUORS
    var REPORT

    var chart
    var dt

    function updateChart(reportData) {
      if (!reportData || !reportData.length) {
        while (chart.series.length > 0) {
          chart.series[0].remove(false)
        }
        chart.redraw()
        return
      }

      // Convert nvd3-style data to Highcharts format
      var series = reportData.map((group, index) => {
        // Get the default Highcharts color for this series
        var seriesColor =
          Highcharts.getOptions().colors[index % Highcharts.getOptions().colors.length]

        return {
          name: group.key,
          data: group.values.map((d) => {
            // Determine symbol shape: circles for proteins, squares for dyes
            var isProtein = FLUORS?.[d.fluor_slug] && FLUORS[d.fluor_slug].type === "p"
            var symbol = isProtein ? "circle" : "square"

            // Scale brightness with a cap to prevent huge symbols
            // Use square root scaling to reduce the range
            var radius = Math.max(2, Math.min(12, Math.sqrt(+d.brightness || 1) * 2))

            // Convert color to rgba with opacity
            var color = Highcharts.color(seriesColor).setOpacity(0.6).get("rgba")

            return {
              x: d.ex_eff && d.ex_eff !== null ? d.ex_eff : 0,
              y: d.em_eff && d.em_eff !== null ? d.em_eff : 0,
              color: color,
              marker: {
                symbol: symbol,
                radius: radius,
              },
              fluor: d.fluor,
              fluor_slug: d.fluor_slug,
              brightness: d.brightness,
              url: d.url,
            }
          }),
        }
      })

      // Remove existing series
      while (chart.series.length > 0) {
        chart.series[0].remove(false)
      }

      // Add new series
      series.forEach((s) => {
        chart.addSeries(s, false)
      })

      chart.redraw()
    }

    function updateData() {
      $.get(`${window.location}json/`, (d) => {
        $("#update-alert").show()
        if (!d.report?.length) {
          $("#status").html('<p class="mt-5">No data yet.  Please update scope report.</p>')
          $("#status").show()
          $(".needs-data").hide()
          $("#report_chart").height(0)
          return
        }

        FLUORS = d.fluors
        REPORT = d.report
        $("#report_chart").height(750)

        updateChart(REPORT)

        $("#status").hide()
        $(".needs-data").show()

        var colHeaders = ["ex_eff", "ex_eff_broad", "em_eff", "brightness"]
        var colTitles = [
          { title: "Fluor", class: "fluor", width: 180, targets: 0 },
          {
            title: "EC",
            class: "col_ext_coeff",
            visible: localStorage.getItem(`${SCOPE_ID}ext_coeff_checkbox`) === "true",
          },
          {
            title: "QY",
            class: "col_qy",
            visible: localStorage.getItem(`${SCOPE_ID}qy_checkbox`) === "true",
          },
          { title: "Ex Max", class: "col_ex_max", visible: false },
          { title: "Em Max", class: "col_em_max", visible: false },
          { title: "Agg", class: "col_agg", visible: false },
          { title: "Switch Type", class: "col_switch_type", visible: false },
          { title: "UUID", class: "col_uuid", visible: false },
        ]
        var ocNames = []
        var fluorData = {}

        d.report.sort(compare)
        // d.report is an array, where each item is an object corresponding to a different OC
        for (let i = d.report.length - 1; i >= 0; i--) {
          const values = d.report[i].values // the array of values for this OC
          const thisOC = d.report[i].key // the name of this thisOC
          ocNames.push(thisOC)

          for (let x = 0; x < colHeaders.length; x++) {
            const cls = `${slugify(`col_${thisOC}`)} col_${colHeaders[x]}`
            const vis =
              localStorage.getItem(`${SCOPE_ID + slugify(thisOC)}_checkbox`) !== "false" &&
              $(`#${colHeaders[x]}_checkbox`).prop("checked")
            colTitles.push({
              title: `${thisOC} ${titleCase(colHeaders[x])}`,
              class: cls,
              visible: vis,
            })
          }

          for (let r = values.length - 1; r >= 0; r--) {
            //extract the values for this OC / Fluor combination
            const thisFluor = values[r].fluor_slug
            if (fluorData[thisFluor] === undefined) {
              fluorData[thisFluor] = { name: values[r].fluor }
            }

            if (fluorData[thisFluor][thisOC] === undefined) {
              fluorData[thisFluor][thisOC] = []
            }

            for (let h = 0; h < colHeaders.length; h++) {
              if (values[r][colHeaders[h]] !== null) {
                fluorData[thisFluor][thisOC].push(values[r][colHeaders[h]])
              } else {
                fluorData[thisFluor][thisOC].push("-")
              }
            }
          }
        }

        // change object into array for Datatables
        var dataRows = []

        var empty = []
        for (let x = colHeaders.length - 1; x >= 0; x--) {
          empty.push("-")
        }

        for (var fluor in fluorData) {
          // skip loop if the property is from prototype
          if (!Object.hasOwn(fluorData, fluor)) continue

          let fluorlink = fluorData[fluor].name
          if (FLUORS[fluor].type === "p") {
            //var fluorlink = '<a href="/protein/'+ FLUORS[fluor].slug.split("_")[0] + '">' + fluorlink + '</a>'
            fluorlink =
              '<a target="_blank" href="' +
              SCOPE_URL +
              "?p=" +
              FLUORS[fluor].slug +
              '">' +
              fluorlink +
              "</a>"
          }
          let f = [
            fluorlink,
            FLUORS[fluor].ext_coeff || "-",
            FLUORS[fluor].qy || "-",
            FLUORS[fluor].ex_max,
            FLUORS[fluor].em_max,
            FLUORS[fluor].agg,
            FLUORS[fluor].switch_type,
            FLUORS[fluor].uuid,
          ]
          for (let o = 0; o < ocNames.length; o++) {
            if (Object.hasOwn(fluorData[fluor], ocNames[o])) {
              f = f.concat(fluorData[fluor][ocNames[o]])
            } else {
              f = f.concat(empty)
            }
          }
          dataRows.push(f)
        }
        $("#oc-toggles").empty()
        for (let o = 0; o < ocNames.length; o++) {
          const ischecked =
            localStorage.getItem(`${SCOPE_ID + slugify(ocNames[o])}_checkbox`) !== "false"
          const el = $("<div>", { class: "custom-control custom-checkbox ml-3" })
            .append(
              $("<input>", {
                class: "custom-control-input toggle-vis",
                type: "checkbox",
                id: `${slugify(ocNames[o])}_checkbox`,
                value: slugify(ocNames[o]),
                checked: ischecked,
              })
            )
            .append(
              $("<label>", {
                class: "custom-control-label",
                for: `${slugify(ocNames[o])}_checkbox`,
              }).text(ocNames[o])
            )

          el.appendTo($("#oc-toggles"))
        }

        var buttonCommon = {
          init: (_api, node, _config) => {
            $(node).removeClass("btn-secondary")
          },
          className: "btn-sm btn-primary",
          exportOptions: {
            format: {
              body: (data, _row, column, node) => {
                switch (column) {
                  case 0:
                    if (data.startsWith("<a")) {
                      return node.childNodes[0].innerHTML
                    }
                    return data
                  default:
                    if (data === "-") {
                      return ""
                    }
                    return data
                }
              },
            },
          },
        }

        colTitles = colTitles.map((d) => {
          d.orderSequence = ["desc", "asc"]
          return d
        })

        dt = $("#report_table").DataTable({
          retrieve: true,
          scrollX: "100%",
          fixedColumns: true,
          columns: colTitles,
          orderSequence: ["desc", "ascending"],
          dom:
            "<'row'<'col-sm-12 col-md-6'l><'col-sm-12 col-md-6'f>>" +
            "<'row'<'col-sm-12'tr>>" +
            "<'row mt-2 small text-small'<'col-sm-12 col-md-3 d-none d-md-block'B><'col-xs-12 col-sm-5 col-md-4'i><'col-xs-12 col-sm-7 col-md-5'p>>",
          buttons: [
            $.extend(true, {}, buttonCommon, {
              extend: "copyHtml5",
            }),
            $.extend(true, {}, buttonCommon, {
              extend: "excelHtml5",
            }),
            $.extend(true, {}, buttonCommon, {
              extend: "csvHtml5",
            }),
            {
              text: "JSON",
              className: "btn-sm btn-primary json-button",
              action: (_e, _dt, _node, _config) => {
                $.fn.dataTable.fileSave(
                  new Blob([JSON.stringify(d.report)]),
                  `report_${SCOPE_ID}.json`
                )
              },
            },
          ],
        })

        dt.clear()
        dt.rows.add(dataRows)
        dt.draw()

        $("#report_table_filter input").on("keyup", function () {
          var rx = this.value.replace(/,\s*$/g, "").replace(/,\s*/g, "|.*")
          dt.search(rx, true).draw()
        })
      })
    }

    let line
    if ($("#update-progress").length) {
      line = new ProgressBar.Line("#update-progress", {
        strokeWidth: 4,
        duration: interval,
        color: "#74779B",
        trailColor: "#ddd",
        trailWidth: 1,
        svgStyle: { width: "100%", height: "100%" },
      })
    }

    function reset_button() {
      clearInterval(window.CHECKER)
      $("#request-report").removeClass("cancel")
      $("#request-report").removeClass("btn-danger")
      $("#request-report").val("Update")
      updateData()
    }

    var LAST_PROG = 0
    var LAST_TIME

    function check_status(next) {
      next = next || interval
      if (JOB_ID !== null) {
        const formData = new URLSearchParams({
          action: "check",
          job_id: JOB_ID,
          csrfmiddlewaretoken: CSRF_TOKEN,
        })

        fetch("", {
          method: "POST",
          headers: {
            "Content-Type": "application/x-www-form-urlencoded",
            // Legacy header required by Django is_ajax() check in dual-purpose endpoints
            "X-Requested-With": "XMLHttpRequest",
          },
          body: formData,
        })
          .then((response) => response.json())
          .then((data) => {
            if (data.ready || data.info === undefined || data.info === null) {
              // job went to completion
              line.animate(1, { duration: 200 }, () => {
                $("#update-alert").fadeOut()
              })
              reset_button()
            } else {
              const passed = Date.now() - LAST_TIME
              const curRate = passed / (data.info.current - LAST_PROG)
              const nextPosition = data.info.current + interval / curRate
              const timeRemaining = (data.info.total - data.info.current) * curRate
              const nextCheck = Math.min(interval, timeRemaining)
              LAST_PROG = data.info.current
              LAST_TIME = Date.now()

              line.animate(Math.min(nextPosition / data.info.total, 1), {
                duration: nextCheck,
              })

              setTimeout(check_status, nextCheck)
              //window.CHECKER = setInterval(check_status, interval);
            }
          })
      }
    }

    $("#request-report").click((e) => {
      $("#request-report").attr("disabled", true)
      $("#request-report").val("Running")

      e.preventDefault()
      //var form = $(this).closest("form")
      let action = "update"
      if ($("#request-report").hasClass("cancel")) {
        action = "cancel"
      }

      const formData = new URLSearchParams({
        action: action,
        scope_id: SCOPE_ID,
        csrfmiddlewaretoken: CSRF_TOKEN,
        job_id: JOB_ID,
        outdated: JSON.stringify(OUTDATED),
      })

      fetch("", {
        method: "POST",
        headers: {
          "Content-Type": "application/x-www-form-urlencoded",
          // Legacy header required by Django is_ajax() check in dual-purpose endpoints
          "X-Requested-With": "XMLHttpRequest",
        },
        body: formData,
      })
        .then((response) => response.json())
        .then((data) => {
          if (data.canceled) {
            line.animate(0, { duration: 200 })
            JOB_ID = null
            reset_button()
          } else if (data.waiting) {
            alert(
              "Too many jobs currently running on the server, please wait a minute and try again."
            )
            JOB_ID = null
            reset_button()
          } else {
            // start request
            line.animate(0, { duration: 10 })
            $("#request-report").addClass("cancel")
            $("#request-report").addClass("btn-danger")
            JOB_ID = data.job
            LAST_TIME = Date.now()
            setTimeout(check_status, 250, interval - 250)
          }
        })
        .catch((_error) => {
          $("#alert-msg").text(
            "Sorry!  There was an unexpected error on the server.  Please try again in a few minutes, or contact us if the problem persists."
          )
          $("#update-alert").removeClass("alert-info")
          $("#update-alert").addClass("alert-danger")
          $("#request-report").val("Error!")
          $("#request-report").addClass("btn-danger")
          setTimeout(() => {
            reset_button()
            $("#request-report").attr("disabled", false)
            $("#alert-msg").text("try again?...")
            $("#update-alert").removeClass("alert-danger")
            $("#update-alert").addClass("alert-info")
            $("#request-report").removeClass("btn-danger")
          }, 20000)
        })
    })

    $(document).ready(() => {
      // create the chart with Highcharts
      chart = Highcharts.chart("report_chart", {
        chart: {
          type: "scatter",
          zoomType: "xy",
        },
        title: { text: null },
        legend: {
          enabled: true,
          align: "center",
          verticalAlign: "top",
          layout: "horizontal",
          itemStyle: {
            fontSize: "12px",
          },
        },
        xAxis: {
          title: { text: "Excitation Efficiency" },
          labels: {
            format: "{value:.2f}",
          },
          min: -0.02,
          max: 1.0,
          gridLineWidth: 1,
          gridLineColor: "#e6e6e6",
          lineWidth: 0,
          tickColor: "#e6e6e6",
        },
        yAxis: {
          title: { text: "Collection Efficiency" },
          labels: {
            format: "{value:.2f}",
          },
          gridLineWidth: 1,
          gridLineColor: "#e6e6e6",
        },
        tooltip: {
          useHTML: true,
          formatter: function () {
            const point = this.point
            const fluor = FLUORS?.[point.fluor_slug]

            const formatPercent = (value) => `${Math.round(value * 1000) / 10}%`
            const formatValue = (value) => value ?? "-"

            const fields = [
              { label: "Ex Eff", value: formatPercent(point.x) },
              { label: "Em Eff", value: formatPercent(point.y) },
              { label: "Brightness", value: formatValue(point.brightness) },
              { label: "EC", value: formatValue(fluor?.ext_coeff) },
              { label: "QY", value: formatValue(fluor?.qy) },
            ]

            const rows = fields.map((f) => `<span>${f.label}: ${f.value}</span><br>`).join("")

            return `<div><p><strong>${point.fluor}</strong></p>${rows}</div>`
          },
        },
        plotOptions: {
          legend: {
            squareSymbol: true,
            symbolRadius: 12,
            symbolHeight: 12,
            symbolWidth: 12,
          },
          series: {
            events: {
              legendItemClick: function () {
                var chart = this.chart
                var currentTime = Date.now()

                // Check for double-click (within 300ms)
                if (this._lastClickTime && currentTime - this._lastClickTime < 300) {
                  // Double-click: show only this series
                  let anyOthersVisible = false
                  chart.series.forEach((series) => {
                    if (series !== this && series.visible) {
                      anyOthersVisible = true
                    }
                  })

                  if (anyOthersVisible) {
                    // Hide all others, show only this one
                    chart.series.forEach((series) => {
                      series.setVisible(series === this, false)
                    })
                  } else {
                    // All others are hidden, restore all
                    chart.series.forEach((series) => {
                      series.setVisible(true, false)
                    })
                  }
                  chart.redraw()
                  this._lastClickTime = null
                  return false // Prevent default toggle
                } else {
                  // Single-click: allow default toggle behavior
                  this._lastClickTime = currentTime
                  return true // Allow default toggle
                }
              },
            },
          },
          scatter: {
            legendSymbol: "rectangle",
            marker: {
              radius: 5,
              lineWidth: 0.5,
              lineColor: "#0000004d",
              states: {
                hover: {
                  enabled: true,
                },
              },
            },
            states: {
              hover: {
                marker: {
                  enabled: false,
                },
              },
            },
            point: {
              events: {
                click: function () {
                  if (this.url) {
                    window.open(this.url)
                  }
                },
              },
            },
          },
        },
        series: [],
        lang: {
          noData: "",
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

      updateData()

      $("body").on("click", "input.toggle-vis, input.meas-vis", function (_e) {
        localStorage.setItem(SCOPE_ID + $(this).attr("id"), $(this).prop("checked"))
        if (["ext_coeff", "qy"].includes($(this).val())) {
          dt.columns(`.col_${$(this).val()}`).visible($(this).prop("checked"))
          return
        }
        if ($(this).prop("checked")) {
          // prevent "re-enabling disabled things"
          let columns = dt.columns(`.col_${$(this).val()}`)[0]
          let included
          if ($(this).hasClass("toggle-vis")) {
            included = $("input.meas-vis:checkbox:checked")
              .map((_x, d) => `.col_${d.value}`)
              .toArray()
              .join(", ")
          } else {
            included = $("input.toggle-vis:checkbox:checked")
              .map((_x, d) => `.col_${d.value}`)
              .toArray()
              .join(", ")
          }
          included = dt.columns(included)[0]
          columns = columns.filter((elem) => included.includes(elem))
          dt.columns(columns).visible($(this).prop("checked"))
        } else {
          dt.columns(`.col_${$(this).val()}`).visible($(this).prop("checked"))
        }
      })

      $(".table-filter").change(function () {
        var uuids = $("#probe_filter option:selected").data("ids")

        // filter the table
        localStorage.setItem(SCOPE_ID + $(this).attr("id"), this.value)
        var searchval = this.value
        if (searchval !== "") {
          searchval = `^${this.value}$`
        }
        if (this.id === "probe_filter" && uuids) {
          searchval = `^(${uuids.join("|")})$`
        }
        var searchcol = $(this).data("col")
        dt.column(`.col_${searchcol}`).search(searchval, true, false).draw()

        // filter the chart
        var newReport = []
        const filterfunc = (d) =>
          ($("#probe_filter").val() === "" ||
            (FLUORS[d.fluor_slug].type === "p" && uuids.includes(FLUORS[d.fluor_slug].uuid))) &&
          ($("#agg_filter").val() === "" || FLUORS[d.fluor_slug].agg === $("#agg_filter").val()) &&
          ($("#switch_filter").val() === "" ||
            FLUORS[d.fluor_slug].switch_type === $("#switch_filter").val())

        for (let i = 0; i < REPORT.length; i++) {
          newReport.push(Object.assign({}, REPORT[i]))
          newReport[i].values = newReport[i].values.filter(filterfunc)
        }

        updateChart(newReport)
      })

      $("#meas-toggles").empty()
      var measbuttons = [
        ["ext_coeff", "EC"],
        ["qy", "QY"],
        ["ex_eff", "Ex Eff"],
        ["ex_eff_broad", "Ex Eff (Broadband)"],
        ["em_eff", "Em Eff"],
        ["brightness", "Brightness"],
      ]
      for (let o = 0; o < measbuttons.length; o++) {
        const id = `${measbuttons[o][0]}_checkbox`
        let ischecked = localStorage.getItem(SCOPE_ID + id)
        if (
          ischecked === null &&
          ["ex_eff_checkbox", "em_eff_checkbox", "brightness_checkbox"].includes(id)
        ) {
          ischecked = true
        } else {
          ischecked = ischecked === "true"
        }

        const el = $("<div>", { class: "custom-control custom-checkbox ml-3" })
          .append(
            $("<input>", {
              class: "custom-control-input meas-vis",
              type: "checkbox",
              id: id,
              value: measbuttons[o][0],
              checked: ischecked,
            })
          )
          .append($("<label>", { class: "custom-control-label", for: id }).text(measbuttons[o][1]))
        el.appendTo($("#meas-toggles"))
      }
    })
  },
})

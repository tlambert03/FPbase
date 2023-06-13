import ProgressBar from "progressbar.js"
import $ from "jquery"
import nv from "nvd3"
import "../css/nv.d3.css"
import d3 from "d3"

function compare(a, b) {
  if (a.key < b.key) return 1
  if (a.key > b.key) return -1
  return 0
}

function titleCase(str) {
  str = str
    .replace(/_/g, " ")
    .toLowerCase()
    .split(" ")
  for (var i = 0; i < str.length; i++) {
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
  scopeReport: function(config) {
    var CSRF_TOKEN = config.csrfToken
    var SCOPE_ID = config.scopeID
    var SCOPE_URL = config.scopeURL
    var JOB_ID = null
    var OUTDATED = undefined
    if (config.outdated.length) {
      OUTDATED = config.outdated
    }
    var interval = config.interval || 3500
    var FLUORS
    var REPORT

    var chart
    var chartData
    var dt

    function updateData() {
      $.get(window.location + "json/", function(d) {
        $("#update-alert").show()
        if (!(d.report && d.report.length)) {
          $("#status").html(
            '<p class="mt-5">No data yet.  Please update scope report.</p>'
          )
          $("#status").show()
          $(".needs-data").hide()
          $("#report_chart").height(0)
          return
        }

        FLUORS = d.fluors
        REPORT = d.report
        $("#report_chart").height(750)
        chartData
          .datum(REPORT)
          .transition()
          .duration(300)
          .call(chart)
        $("#status").hide()
        $(".needs-data").show()

        var colHeaders = ["ex_eff", "ex_eff_broad", "em_eff", "brightness"]
        var colTitles = [
          { title: "Fluor", class: "fluor", width: 180, targets: 0 },
          {
            title: "EC",
            class: "col_ext_coeff",
            visible:
              localStorage.getItem(SCOPE_ID + "ext_coeff_checkbox") === "true"
          },
          {
            title: "QY",
            class: "col_qy",
            visible: localStorage.getItem(SCOPE_ID + "qy_checkbox") === "true"
          },
          { title: "Ex Max", class: "col_ex_max", visible: false },
          { title: "Em Max", class: "col_em_max", visible: false },
          { title: "Agg", class: "col_agg", visible: false },
          { title: "Switch Type", class: "col_switch_type", visible: false },
          { title: "UUID", class: "col_uuid", visible: false }
        ]
        var ocNames = []
        var fluorData = {}

        d.report.sort(compare)
        // d.report is an array, where each item is an object corresponding to a different OC
        for (var i = d.report.length - 1; i >= 0; i--) {
          var values = d.report[i].values // the array of values for this OC
          var thisOC = d.report[i].key // the name of this thisOC
          ocNames.push(thisOC)

          for (var x = 0; x < colHeaders.length; x++) {
            var cls = slugify("col_" + thisOC) + " col_" + colHeaders[x]
            var vis =
              localStorage.getItem(SCOPE_ID + slugify(thisOC) + "_checkbox") !==
                "false" && $("#" + colHeaders[x] + "_checkbox").prop("checked")
            colTitles.push({
              title: thisOC + " " + titleCase(colHeaders[x]),
              class: cls,
              visible: vis
            })
          }

          for (var r = values.length - 1; r >= 0; r--) {
            //extract the values for this OC / Fluor combination
            var thisFluor = values[r].fluor_slug
            if (fluorData[thisFluor] === undefined) {
              fluorData[thisFluor] = { name: values[r].fluor }
            }

            if (fluorData[thisFluor][thisOC] === undefined) {
              fluorData[thisFluor][thisOC] = []
            }

            for (var h = 0; h < colHeaders.length; h++) {
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
          if (!fluorData.hasOwnProperty(fluor)) continue

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
          var f = [
            fluorlink,
            FLUORS[fluor]["ext_coeff"] || "-",
            FLUORS[fluor]["qy"] || "-",
            FLUORS[fluor]["ex_max"],
            FLUORS[fluor]["em_max"],
            FLUORS[fluor]["agg"],
            FLUORS[fluor]["switch_type"],
            FLUORS[fluor]["uuid"]
          ]
          for (var o = 0; o < ocNames.length; o++) {
            if (fluorData[fluor].hasOwnProperty(ocNames[o])) {
              f = f.concat(fluorData[fluor][ocNames[o]])
            } else {
              f = f.concat(empty)
            }
          }
          dataRows.push(f)
        }
        $("#oc-toggles").empty()
        for (let o = 0; o < ocNames.length; o++) {
          var ischecked =
            localStorage.getItem(
              SCOPE_ID + slugify(ocNames[o]) + "_checkbox"
            ) !== "false"
          var el = $("<div>", { class: "custom-control custom-checkbox ml-3" })
            .append(
              $("<input>", {
                class: "custom-control-input toggle-vis",
                type: "checkbox",
                id: slugify(ocNames[o]) + "_checkbox",
                value: slugify(ocNames[o]),
                checked: ischecked
              })
            )
            .append(
              $("<label>", {
                class: "custom-control-label",
                for: slugify(ocNames[o]) + "_checkbox"
              }).text(ocNames[o])
            )

          el.appendTo($("#oc-toggles"))
        }

        var buttonCommon = {
          init: function(api, node, config) {
            $(node).removeClass("btn-secondary")
          },
          className: "btn-sm btn-primary",
          exportOptions: {
            format: {
              body: function(data, row, column, node) {
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
              }
            }
          }
        }

        colTitles = colTitles.map(function(d) {
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
              extend: "copyHtml5"
            }),
            $.extend(true, {}, buttonCommon, {
              extend: "excelHtml5"
            }),
            $.extend(true, {}, buttonCommon, {
              extend: "csvHtml5"
            }),
            {
              text: "JSON",
              className: "btn-sm btn-primary json-button",
              action: function(e, dt, node, config) {
                $.fn.dataTable.fileSave(
                  new Blob([JSON.stringify(d.report)]),
                  "report_" + SCOPE_ID + ".json"
                )
              }
            }
          ]
        })

        dt.clear()
        dt.rows.add(dataRows)
        dt.draw()

        $("#report_table_filter input").on("keyup", function() {
          var rx = this.value.replace(/,\s*$/g, "").replace(/,\s*/g, "|.*")
          try {
            //var regex = new RegExp(rx)
          } catch (e) {
            console.log(e)
            return
            // don't need to know if this regex fails...
          }
          dt.search(rx, true).draw()
        })
      })
    }

    if ($("#update-progress").length) {
      var line = new ProgressBar.Line("#update-progress", {
        strokeWidth: 4,
        duration: interval,
        color: "#74779B",
        trailColor: "#ddd",
        trailWidth: 1,
        svgStyle: { width: "100%", height: "100%" }
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
        $.post({
          url: "",
          data: {
            action: "check",
            job_id: JOB_ID,
            csrfmiddlewaretoken: CSRF_TOKEN
          },
          dataType: "json",
          success: function(data) {
            if (data.ready || data.info === undefined || data.info === null) {
              // job went to completion
              line.animate(1, { duration: 200 }, function() {
                $("#update-alert").fadeOut()
              })
              reset_button()
            } else {
              var passed = Date.now() - LAST_TIME
              var curRate = passed / (data.info.current - LAST_PROG)
              var nextPosition = data.info.current + interval / curRate
              var timeRemaining =
                (data.info.total - data.info.current) * curRate
              var nextCheck = Math.min(interval, timeRemaining)
              LAST_PROG = data.info.current
              LAST_TIME = Date.now()

              line.animate(Math.min(nextPosition / data.info.total, 1), {
                duration: nextCheck
              })

              setTimeout(check_status, nextCheck)
              //window.CHECKER = setInterval(check_status, interval);
            }
          }
        })
      }
    }

    $("#request-report").click(function(e) {
      $("#request-report").attr("disabled", true)
      $("#request-report").val("Running")

      e.preventDefault()
      //var form = $(this).closest("form")
      let action = "update"
      if ($("#request-report").hasClass("cancel")) {
        action = "cancel"
      }
      $.post({
        url: "",
        data: {
          action: action,
          scope_id: SCOPE_ID,
          csrfmiddlewaretoken: CSRF_TOKEN,
          job_id: JOB_ID,
          outdated: JSON.stringify(OUTDATED)
        },
        dataType: "json",
        success: function(data) {
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
        },
        error: function(req, status, err) {
          $("#alert-msg").text(
            "Sorry!  There was an unexpected error on the server.  Please try again in a few minutes, or contact us if the problem persists."
          )
          $("#update-alert").removeClass("alert-info")
          $("#update-alert").addClass("alert-danger")
          $("#request-report").val("Error!")
          $("#request-report").addClass("btn-danger")
          setTimeout(function() {
            reset_button()
            $("#request-report").attr("disabled", false)
            $("#alert-msg").text("try again?...")
            $("#update-alert").removeClass("alert-danger")
            $("#update-alert").addClass("alert-info")
            $("#request-report").removeClass("btn-danger")
          }, 20000)
        }
      })
    })

    $(document).ready(function() {
      // create the chart
      nv.addGraph(function() {
        chart = nv.models
          .scatterChart()
          .x(function(d) {
            if (d.ex_eff && d.ex_eff !== null) {
              return d.ex_eff
            }
            return 0
          })
          .y(function(d) {
            if (d.em_eff && d.em_eff !== null) {
              return d.em_eff
            }
            return 0
          })
          .pointSize(function(d) {
            return Math.max(1, +d.brightness * 2)
          })
          .pointRange([10, 400])
          .pointDomain([1, 80])
          //.showDistX(true)
          //.showDistY(true)
          //.yDomain([0, 1])
          //.xDomain([0, 1])
          .useVoronoi(true)
          .noData("")

        chart.xAxis
          .tickFormat(d3.format(".02f"))
          .axisLabel("Excitation Efficiency")
        chart.yAxis
          .tickFormat(d3.format(".02f"))
          .axisLabel("Collection Efficiency")

        chartData = d3.select("#report_chart svg")
        chartData
          .datum([])
          .transition()
          .duration(300)
          .call(chart)
        nv.utils.windowResize(chart.update)
        chart.tooltip.contentGenerator(function(key, x, y, e, graph) {
          return (
            '<div class="report-tooltip" style="margin-top:' +
            window.scrollY +
            'px;"><p><strong>' +
            key.point.fluor +
            "</strong></p>" +
            "<p>Ex Eff: " +
            Math.round(1000 * key.point.x) / 10 +
            " %</p>" +
            "<p>Em Eff: " +
            Math.round(1000 * key.point.y) / 10 +
            " %</p>" +
            "<p>Brightness: " +
            (key.point.brightness || "-") +
            "</p>" +
            "<p>EC: " +
            (FLUORS[key.point.fluor_slug]["ext_coeff"] || "-") +
            "</p>" +
            "<p>QY: " +
            (FLUORS[key.point.fluor_slug]["qy"] || "-") +
            "</p>" +
            "</div>"
          )
        })
        chart.scatter.dispatch.on("elementClick", function(event) {
          //location.href = event.point.url;
          window.open(event.point.url)
        })
        return chart
      })
      updateData()

      $("body").on("click", "input.toggle-vis, input.meas-vis", function(e) {
        localStorage.setItem(
          SCOPE_ID + $(this).attr("id"),
          $(this).prop("checked")
        )
        if (["ext_coeff", "qy"].includes($(this).val())) {
          dt.columns(".col_" + $(this).val()).visible($(this).prop("checked"))
          return
        }
        if ($(this).prop("checked")) {
          // prevent "re-enabling disabled things"
          var columns = dt.columns(".col_" + $(this).val())[0]
          let included
          if ($(this).hasClass("toggle-vis")) {
            included = $("input.meas-vis:checkbox:checked")
              .map(function(x, d) {
                return ".col_" + d.value
              })
              .toArray()
              .join(", ")
          } else {
            included = $("input.toggle-vis:checkbox:checked")
              .map(function(x, d) {
                return ".col_" + d.value
              })
              .toArray()
              .join(", ")
          }
          included = dt.columns(included)[0]
          columns = columns.filter(function(elem) {
            return included.includes(elem)
          })
          dt.columns(columns).visible($(this).prop("checked"))
        } else {
          dt.columns(".col_" + $(this).val()).visible($(this).prop("checked"))
        }
      })

      $(".table-filter").change(function() {
        var uuids = $("#probe_filter option:selected").data("ids")

        // filter the table
        localStorage.setItem(SCOPE_ID + $(this).attr("id"), this.value)
        var searchval = this.value
        if (searchval !== "") {
          searchval = "^" + this.value + "$"
        }
        if (this.id === "probe_filter" && uuids) {
          searchval = "^(" + uuids.join("|") + ")$"
        }
        var searchcol = $(this).data("col")
        dt.column(".col_" + searchcol)
          .search(searchval, true, false)
          .draw()

        // filter the chart
        var newReport = []
        const filterfunc = d =>
          ($("#probe_filter").val() === "" ||
            (FLUORS[d.fluor_slug].type === "p" &&
              uuids.includes(FLUORS[d.fluor_slug].uuid))) &&
          ($("#agg_filter").val() === "" ||
            FLUORS[d.fluor_slug].agg === $("#agg_filter").val()) &&
          ($("#switch_filter").val() === "" ||
            FLUORS[d.fluor_slug].switch_type === $("#switch_filter").val())

        for (var i = 0; i < REPORT.length; i++) {
          newReport.push(Object.assign({}, REPORT[i]))
          newReport[i].values = newReport[i].values.filter(filterfunc)
        }

        d3.select("#report_chart svg").datum(newReport)
        chart.update()
      })

      $("#meas-toggles").empty()
      var measbuttons = [
        ["ext_coeff", "EC"],
        ["qy", "QY"],
        ["ex_eff", "Ex Eff"],
        ["ex_eff_broad", "Ex Eff (Broadband)"],
        ["em_eff", "Em Eff"],
        ["brightness", "Brightness"]
      ]
      for (var o = 0; o < measbuttons.length; o++) {
        var id = measbuttons[o][0] + "_checkbox"
        var ischecked = localStorage.getItem(SCOPE_ID + id)
        if (
          ischecked === null &&
          [
            "ex_eff_checkbox",
            "em_eff_checkbox",
            "brightness_checkbox"
          ].includes(id)
        ) {
          ischecked = true
        } else {
          ischecked = ischecked === "true"
        }

        var el = $("<div>", { class: "custom-control custom-checkbox ml-3" })
          .append(
            $("<input>", {
              class: "custom-control-input meas-vis",
              type: "checkbox",
              id: id,
              value: measbuttons[o][0],
              checked: ischecked
            })
          )
          .append(
            $("<label>", { class: "custom-control-label", for: id }).text(
              measbuttons[o][1]
            )
          )
        el.appendTo($("#meas-toggles"))
      }
    })
  }
})

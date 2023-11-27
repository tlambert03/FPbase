// to anyone reading this... I apologize for how crapy this code is...
// i intend to refactor the whole program... but this was one of those times
// where a bad foundation just kept going.  the program works relatively well,
// but it sure isn't pretty.

import noUiSlider from "nouislider"
import $ from "jquery"
import nv from "nvd3"
import "../css/nv.d3.css"
import d3 from "d3"

;

(function() {
  if (!$(".microscope-wrapper #spectra svg").length) {
    return
  }

  // Safari 3.0+ "[object HTMLElementConstructor]"
  const isSafari =
    /constructor/i.test(window.HTMLElement) ||
    (function(p) {
      return p.toString() === "[object SafariRemoteNotification]"
    })(
      !window.safari ||
        (typeof safari !== "undefined" && window.safari.pushNotification)
    )

  if (!String.prototype.endsWith) {
    String.prototype.endsWith = function(search, this_len) {
      if (typeof this_len === "undefined" || this_len > this.length) {
        this_len = this.length
      }
      return this.substring(this_len - search.length, this_len) === search
    }
  }

  const urlSegments = window.location.pathname.split('/');
  const ScopeID = urlSegments[urlSegments.length - 2]; // The segment before the last '/'
  const KeyPrefix = `microscope_${ScopeID}`;

  let chart
  const CONST = {
    category: {
      dye: "d",
      protein: "p",
      light: "l",
      filter: "f",
      camera: "c",
    },
    stype: {
      ex: "ex",
      abs: "ab",
      em: "em",
      twop: "2p",
      bp: "bp",
      bpx: "bx",
      bpm: "bm",
      sp: "sp",
      lp: "lp",
      bs: "bs",
      qe: "qe",
      pd: "pd",
    },
  }
  let data = []
  const localData = {}
  const options = {
    showArea: localStorage.getItem(`${KeyPrefix}_showArea`) !== "false",
    stickySpectra: localStorage.getItem(`${KeyPrefix}_stickySpectra`),
    minwave: 350,
    maxwave: 800,
    startingBrush: [350, 800],
    // autoscaleBrush: false,
    exNormWave: undefined,
    scale: "linear",
    hide2p: true,
    normMergedEx: localStorage.getItem(`${KeyPrefix}_normMergedEx`) !== "false",
    normMergedScalar: 1,
    focusEnable: localStorage.getItem(`${KeyPrefix}_focusEnable`) !== "false",
    //    scaleToEC: false,
    scaleToQY: localStorage.getItem(`${KeyPrefix}_scaleToQY`) === "true",
    oneAtaTime: localStorage.getItem(`${KeyPrefix}_oneAtaTime`) !== "false",
    precision:
      localStorage.getItem(`${KeyPrefix}_precision`) || (isSafari ? 2 : 1),
    interpolate: localStorage.getItem(`${KeyPrefix}_interpolate`) === "true",
    // modified later at page load
    calcEff: true,
    exEffBroadband:
      localStorage.getItem(`${KeyPrefix}_exEffBroadband`) === "true",
  }
  const userOptions = {
    calcEff: {
      type: "checkbox",
      msg: "Calculate efficiency on update. (may be slower)",
    },
    oneAtaTime: {
      type: "checkbox",
      msg: "Uncheck other similar filters when selecting a filter",
    },
    normMergedEx: {
      type: "checkbox",
      msg: "Normalize merged excitation and light source spectra",
    },
    showArea: { type: "checkbox", msg: "Fill area under curve" },
    // autoscaleBrush: {type: 'checkbox', msg: 'Auto-rescale X-axis (using zoom above auto-disables this)'},
    //    hide2p: {type: 'checkbox', msg: 'Hide 2-photon spectra by default'},
    //    scaleToEC: {type: 'checkbox', msg: 'Scale excitation spectra to extinction coefficient (% of highest fluor)'},
    focusEnable: { type: "checkbox", msg: "Enable zoom/pan bar" },
    interpolate: { type: "checkbox", msg: "Spline interpolation on chart" },
    precision: {
      type: "number",
      msg: "Precision of chart in nm (higher = faster)",
      max: 8,
      min: 1,
    },
    scaleToQY: {
      type: "checkbox",
      msg: "Scale emission spectra to quantum yield",
    },
    exEffBroadband: {
      type: "checkbox",
      msg:
        '"Broadband" excitation efficiency mode <a href="https://help.fpbase.org/tools/microscopes/efficiency#broadband" target="_blank"><i class="fa fa-question-circle text-muted"></i></a>',
    },
  }

  const mobilecheck = function() {
    let check = false
    ;(function(a) {
      if (
        /(android|bb\d+|meego).+mobile|avantgo|bada\/|blackberry|blazer|compal|elaine|fennec|hiptop|iemobile|ip(hone|od)|iris|kindle|lge |maemo|midp|mmp|mobile.+firefox|netfront|opera m(ob|in)i|palm( os)?|phone|p(ixi|re)\/|plucker|pocket|psp|series(4|6)0|symbian|treo|up\.(browser|link)|vodafone|wap|windows ce|xda|xiino/i.test(
          a
        ) ||
        /1207|6310|6590|3gso|4thp|50[1-6]i|770s|802s|a wa|abac|ac(er|oo|s\-)|ai(ko|rn)|al(av|ca|co)|amoi|an(ex|ny|yw)|aptu|ar(ch|go)|as(te|us)|attw|au(di|\-m|r |s )|avan|be(ck|ll|nq)|bi(lb|rd)|bl(ac|az)|br(e|v)w|bumb|bw\-(n|u)|c55\/|capi|ccwa|cdm\-|cell|chtm|cldc|cmd\-|co(mp|nd)|craw|da(it|ll|ng)|dbte|dc\-s|devi|dica|dmob|do(c|p)o|ds(12|\-d)|el(49|ai)|em(l2|ul)|er(ic|k0)|esl8|ez([4-7]0|os|wa|ze)|fetc|fly(\-|_)|g1 u|g560|gene|gf\-5|g\-mo|go(\.w|od)|gr(ad|un)|haie|hcit|hd\-(m|p|t)|hei\-|hi(pt|ta)|hp( i|ip)|hs\-c|ht(c(\-| |_|a|g|p|s|t)|tp)|hu(aw|tc)|i\-(20|go|ma)|i230|iac( |\-|\/)|ibro|idea|ig01|ikom|im1k|inno|ipaq|iris|ja(t|v)a|jbro|jemu|jigs|kddi|keji|kgt( |\/)|klon|kpt |kwc\-|kyo(c|k)|le(no|xi)|lg( g|\/(k|l|u)|50|54|\-[a-w])|libw|lynx|m1\-w|m3ga|m50\/|ma(te|ui|xo)|mc(01|21|ca)|m\-cr|me(rc|ri)|mi(o8|oa|ts)|mmef|mo(01|02|bi|de|do|t(\-| |o|v)|zz)|mt(50|p1|v )|mwbp|mywa|n10[0-2]|n20[2-3]|n30(0|2)|n50(0|2|5)|n7(0(0|1)|10)|ne((c|m)\-|on|tf|wf|wg|wt)|nok(6|i)|nzph|o2im|op(ti|wv)|oran|owg1|p800|pan(a|d|t)|pdxg|pg(13|\-([1-8]|c))|phil|pire|pl(ay|uc)|pn\-2|po(ck|rt|se)|prox|psio|pt\-g|qa\-a|qc(07|12|21|32|60|\-[2-7]|i\-)|qtek|r380|r600|raks|rim9|ro(ve|zo)|s55\/|sa(ge|ma|mm|ms|ny|va)|sc(01|h\-|oo|p\-)|sdk\/|se(c(\-|0|1)|47|mc|nd|ri)|sgh\-|shar|sie(\-|m)|sk\-0|sl(45|id)|sm(al|ar|b3|it|t5)|so(ft|ny)|sp(01|h\-|v\-|v )|sy(01|mb)|t2(18|50)|t6(00|10|18)|ta(gt|lk)|tcl\-|tdg\-|tel(i|m)|tim\-|t\-mo|to(pl|sh)|ts(70|m\-|m3|m5)|tx\-9|up(\.b|g1|si)|utst|v400|v750|veri|vi(rg|te)|vk(40|5[0-3]|\-v)|vm40|voda|vulc|vx(52|53|60|61|70|80|81|83|85|98)|w3c(\-| )|webc|whit|wi(g |nc|nw)|wmlb|wonu|x700|yas\-|your|zeto|zte\-/i.test(
          a.substr(0, 4)
        )
      )
        check = true
    })(navigator.userAgent || navigator.vendor || window.opera)
    return check
  }

  const chartOptions = function() {
    return {
      focusEnable: options.focusEnable,
      focusShowAxisX: false,
      interpolate: options.interpolate ? "basis" : "linear",
      noData: "Add spectra below ...",
      showLegend: true,
      showXAxis: true,
      showYAxis: true,
      duration: 200,
      useInteractiveGuideline: !mobilecheck(),
      clipEdge: true,
      margin: {
        left: 40,
        bottom: options.focusEnable ? 15 : 58,
      },
      yDomain: [0, 1],
      // xDomain: [350, 800],
      // forceY: [0,1.04],
      // forceX: [350, 800],
    }
  }

  const svg = d3.select(".microscope-wrapper #spectra svg")

  function getData(slug) {
    const dfd = $.Deferred()
    // download if not already downloaded
    if (!(slug in localData)) {
      $.getJSON(`/spectra/${slug}`)
        .done(function(d) {
          for (let n = 0; n < d.spectra.length; n++) {
            d.spectra[n] = padDataLimits(d.spectra[n])
            d.spectra[n].exNormed = 1
            if ((d.spectra[n].type == CONST.stype.twop) & options.hide2p) {
              d.spectra[n].disabled = true
            } else if (
              (d.spectra[n].type == CONST.stype.twop) &
              !options.hide2p
            ) {
              d.spectra[n].disabled = false
            }
            if (
              (d.spectra[n].category == CONST.category.light) |
              (d.spectra[n].category == CONST.category.camera)
            ) {
              d.spectra[n].gradient = true
            }
          }
          localData[slug] = d.spectra
          dfd.resolve(localData[slug])
        })
        .fail(function(d) {
          dfd.reject(d.status)
        })
    } else {
      // otherwise pull from local dict
      dfd.resolve(localData[slug])
    }
    return dfd.promise()
  }

  function padDataLimits(d) {
    for (let i = d.minwave - 1; i >= Math.max(options.minwave, 350); i--) {
      d.values.unshift({ x: i, y: 0 })
    }
    for (let n = d.maxwave + 1; n <= Math.min(options.maxwave, 800); n++) {
      d.values.push({ x: n, y: 0 })
    }
    return d
  }

  function dataHasKey(key) {
    return (
      $.grep(data, function(obj) {
        return obj.key == key
      }).length > 0
    )
  }

  function dataHasSlug(slug) {
    return (
      $.grep(data, function(obj) {
        return obj.slug == slug
      }).length > 0
    )
  }

  function dataItemMatching(filter, d) {
    d = d || data
    return d.filter(function(item) {
      for (const key in filter) {
        if (typeof item[key] === "undefined" || item[key] != filter[key])
          return false
      }
      return true
    })
  }

  function asyncFunction(item) {
    itemsProcessed++
    if (itemsProcessed === array.length) {
      callback()
    }
  }

  function mergeSpectra(concatslug) {
    const dfd = $.Deferred()
    // download if not already downloaded
    if (!(concatslug in localData)) {
      const sp = []
      const slugs = concatslug.split("/")
      slugs.shift() // get rid of "merged"
      if (concatslug.endsWith("normed")) slugs.pop()
      const asyncFunction = function(s, onfinished) {
        if (s.startsWith("laser")) {
          if (localData.hasOwnProperty(s)) {
            sp.push.apply(sp, localData[s])
          } else {
            const F = customLaser(+s.replace("laser-", ""))
            localData[s] = [F]
            sp.push.apply(sp, localData[s])
          }
          onfinished()
        } else if (!s.startsWith("merged")) {
          if (localData.hasOwnProperty(s)) {
            sp.push.apply(sp, localData[s])
            onfinished()
          } else {
            getData(s).then(function() {
              sp.push.apply(sp, localData[s])
              onfinished()
            })
          }
        }
      }
      let itemsProcessed = 0
      slugs.forEach(function(s, ind, array) {
        asyncFunction(s, function() {
          itemsProcessed += 1
          if (itemsProcessed === array.length) {
            const names = slugs
              .map(function(o) {
                return localData[o][0].key
              })
              .reduce(function(a, b) {
                return `${a}+${b}`
              })
            let combined = combineSpectra(sp)
            let _slug = concatslug
            let normmax = 1
            if (options.normMergedEx) {
              const temp = normalizeSpectrum(combined)
              combined = temp[0]
              normmax = temp[1]
              options.normMergedScalar = normmax
              _slug += "/normed"
            }
            localData[concatslug] = [
              {
                slug: _slug,
                key: names,
                values: combined,
                area: true,
                scalar: normmax,
                color: "url(#wavecolor_gradient)",
                classed: "combined-light-ex",
              },
            ]
            dfd.resolve(localData[concatslug])
          }
        })
      })
    } else {
      // otherwise pull from local dict
      options.normMergedScalar = localData[concatslug][0].scalar
      dfd.resolve(localData[concatslug])
    }
    return dfd.promise()
  }

  function pushData(item) {
    if (options.precision > 1) {
      const newitem = { ...item }
      newitem.values = newitem.values.filter(function(d) {
        return d.x % options.precision === 0
      })
      data.push(newitem)
    } else {
      data.push(item)
    }
  }

  function addMerged(comboslug) {
    return mergeSpectra(comboslug).then(function(d) {
      pushData(d[0])
    })
  }

  function addItem(slug, subtype) {
    // add spectral data to array
    subtype = subtype || false

    return getData(slug)
      .then(function(d) {
        for (let i = 0; i < d.length; i++) {
          if ((!subtype | (d[i].type == subtype)) & !dataHasKey(d[i].key)) {
            pushData(JSON.parse(JSON.stringify(d[i]))) // make a copy of the object
          }
        }
      })
      .fail(function(d) {
        console.log(`item not found: ${slug}`)
      })
  }

  function removeSubtypes(subtype) {
    const ix = $.map(data, function(obj, index) {
      if (obj.type == subtype) {
        return index
      }
    })
    for (let i = ix.length - 1; i >= 0; i--) data.splice(ix[i], 1)
  }

  function removeItem(slug, subtype) {
    // optional subtype argument can be used to remove a specific type of spectrum
    // rather than the whole slug family
    subtype = subtype || false

    const ix = $.map(data, function(obj, index) {
      if (obj.slug == slug) {
        if (!subtype | (obj.type == subtype)) {
          return index
        }
      }
    })
    for (let i = ix.length - 1; i >= 0; i--) {
      data.splice(ix[i], 1)
    }
  }

  // /// chart setup

  function refreshChart() {
    chart.lines.duration(200)
    // if (options.autoscaleBrush) {
    //     var smin = 10000;
    //     var smax = 0;
    //     for (var i = 0; i < data.length; i++) {
    //         if (!data[i].disabled) {
    //             smin = Math.min(data[i].minwave, smin);
    //             smax = Math.max(data[i].maxwave, smax);
    //         }
    //     }
    //     chart.brushExtent([smin, smax]);
    // }

    scaleDataToOptions()
    excitationNormalization()
    chart.update()
    updateGlobalGradient()
    chart.lines.duration(0)
    if (options.scaleToQY || options.scaleToEC || options.exNormWave) {
      $("#y-zoom-slider").show()
      resizeYSlider()
    } else {
      $("#y-zoom-slider").hide()
    }
    if ($(document).width() < 576) {
      $("svg.nvd3-svg").height(chart.legend.height() + 210)
      chart.update()
    }
  }

  // // Scaling functions

  function unscale_all() {
    options.exNormWave = undefined
    options.scaleToQY = false
    options.scaleToEC = false
    $("#exnormRadioOFF").prop("checked", true)
    $("#scaleToQY_input").prop("checked", false)
    $("#scaleToEC_input").prop("checked", false)
    refreshChart()
  }

  function scale_data_up(filter) {
    // scale data "up" according to the data.scalar value
    // filter can be .e.g. {type: 'ex'} to scale by ExtCoeff
    const maxScalar =
      Math.max.apply(
        null,
        data.map(function(e) {
          return e.scalar || 0
        })
      ) || 1
    for (let n = 0; n < data.length; n++) {
      // only scale certain data by filter
      let skip = false
      if (data[n].scaled || typeof data[n].scalar === "undefined") {
        skip = true
      }
      for (const key in filter) {
        if (
          typeof data[n][key] === "undefined" ||
          data[n][key] != filter[key]
        ) {
          skip = true
          break
        }
      }

      if (!skip) {
        let SCALE = data[n].scalar || 0.001
        if (data[n].type == "ex") {
          SCALE /= maxScalar
        }
        // do the scaling
        for (let i = 0; i < data[n].values.length; i++) {
          data[n].values[i].y = data[n].values[i].y * SCALE
        }
        data[n].scaled = SCALE
      }
    }
  }

  function scale_data_down(filter) {
    for (let n = 0; n < data.length; n++) {
      // only scale certain data by filter
      let skip = false
      for (const key in filter) {
        if (
          typeof data[n][key] === "undefined" ||
          data[n][key] != filter[key] ||
          !data[n].scaled
        ) {
          skip = true
          break
        }
      }
      if (!skip) {
        for (let i = 0; i < data[n].values.length; i++) {
          data[n].values[i].y = data[n].values[i].y / data[n].scaled
        }
        data[n].scaled = false
      }
    }
  }

  function excitationNormalization() {
    if (options.exNormWave) {
      for (var n = 0; n < data.length; n++) {
        if (data[n].type == "em") {
          const exvals = dataItemMatching({ slug: data[n].slug, type: "ex" })[0]
            .values
          let targetScalar = dataItemMatching(
            { x: options.exNormWave },
            exvals
          )[0].y
          targetScalar = Math.max(targetScalar, 0.0001) // so we don't clobber the values with zero
          const currentScalar = data[n].exNormed
          if (currentScalar != targetScalar) {
            var S = targetScalar / currentScalar
            data[n].values = data[n].values.map(function(item) {
              return { x: item.x, y: item.y * S }
            })
            data[n].exNormed = targetScalar
          }
        }
      }
    } else {
      for (var n = 0; n < data.length; n++) {
        if (data[n].type == "em" && data[n].exNormed != 1) {
          data[n].values = data[n].values.map(function(item) {
            return { x: item.x, y: item.y / data[n].exNormed }
          })
          data[n].exNormed = 1
        }
      }
    }
  }

  function resizeYSlider() {
    try {
      const rect = $("#spectrasvg")
        .find(".nv-focus")
        .first()
        .find(".nv-background>rect")
      const h = +rect.attr("height")
      // var t = Math.abs(rect.offset().top - $('#spectrasvg').offset().top);
      $("#y-zoom-slider").height(h * 0.88)
    } catch (e) {}
  }

  function scaleDataToOptions() {
    if (options.scaleToEC) {
      scale_data_up({ type: "ex" })
    } else {
      scale_data_down({ type: "ex" })
    }

    if (options.scaleToQY) {
      scale_data_up({ type: "em" })
    } else {
      scale_data_down({ type: "em" })
    }
  }

  // // ON WINDOW LOAD
  function getUrlParams(prop) {
    const params = {}
    const search = decodeURIComponent(
      window.location.href.slice(window.location.href.indexOf("?") + 1)
    )
    const definitions = search.split("&")

    definitions.forEach(function(val, key) {
      const parts = val.split("=", 2)
      params[parts[0]] = parts[1]
    })

    return prop && prop in params ? params[prop] : params
  }

  $(".2pcheck").change(function(e) {
    // when the ex/em/2p checkbox is checked,
    // enable/disable the corresponding item in the data
    const slug = $(this)
      .closest(".row")
      .find("select")
      .val()
    const type = $(this).val()
    for (let i = 0; i < data.length; i++) {
      if (data[i].slug == slug && data[i].type == type) {
        data[i].disabled = !this.checked
      }
      if (type == CONST.stype.twop && localData[slug]) {
        const pos = localData[slug]
          .map(function(e) {
            return e.type
          })
          .indexOf("2p")
        localData[slug][pos].disabled = !this.checked
      }
    }
    refreshChart()
  })

  function setup_from_url(urlParams) {
    if ("l" in urlParams) {
      $("#light-select").val(urlParams.l)
    }
    if ("d" in urlParams) {
      $("#camera-select").val(urlParams.d)
    }
    if ("p" in urlParams) {
      // $("#fluor-select").val(urlParams.p).trigger('change.select2');
      $("#fluor-select")
        .val(urlParams.p)
        .change()
    }
    if ("c" in urlParams) {
      $("#config-select option").each(function() {
        if (this.value == urlParams.c) {
          $(this)
            .val(urlParams.c)
            .prop("selected", true)
            .change()
        }
      })
    } else if ("f" in urlParams) {
      $.each(urlParams.f.split(","), function(index, val) {
        if (val.endsWith("!")) {
          val = val.slice(0, -1)
          $(`#${val}`)
            .closest("li")
            .find(".invswitch")
            .prop("checked", true)
        }
        $(`#${val}`).prop("checked", true)
      })
    } else {
      $("#config-select option")
        .eq(1)
        .prop("selected", true)
        .change()
    }

    if ("sticky" in urlParams) {
      options.stickySpectra = Boolean(
        urlParams.sticky === "true" || +urlParams.sticky > 0
      )
    }
    if (options.stickySpectra === true) {
      $("#stickyPin").click()
    }

    setTimeout(updateChart, 300)
  }

  // on page load, setup the chart
  $(function() {
    $("#y-zoom-slider").hide()
    // $('[data-toggle="popover"]').popover()

    const urlParams = getUrlParams()
    options.precision = urlParams.precision || options.precision

    // Set options first from database (populated in microscope_detail.html)
    if (typeof scopecfg !== "undefined") {
      options.calcEff = scopecfg.calcEff;
      options.minwave = scopecfg.minwave;
      options.maxwave = scopecfg.maxwave;
      options.focusEnable = scopecfg.focusEnable;
      options.showArea = scopecfg.showArea;
    }
    // But allow localStorage to override them

    if (`${KeyPrefix}_calcEff` in localStorage) {
      options.calcEff = localStorage.getItem(`${KeyPrefix}_calcEff`) !== "false";
    }
    if (`${KeyPrefix}_minwave` in localStorage) {
      options.minwave = localStorage.getItem(`${KeyPrefix}_minwave`);
    }
    if (`${KeyPrefix}_maxwave` in localStorage) {
      options.maxwave = localStorage.getItem(`${KeyPrefix}_maxwave`);
    }
    if (`${KeyPrefix}_focusEnable` in localStorage) {
      options.focusEnable = localStorage.getItem(`${KeyPrefix}_focusEnable`) !== "false";
    }

    // Then, finally, allow URL params to override them
    if (urlParams.eff !== undefined) {
      options.calcEff = !(urlParams.eff == "false")
    }

    $.each(userOptions, function(key, value) {
      let divClasses =
        value.type === "checkbox" ? "custom-control custom-checkbox" : ""
      divClasses += " mb-1 pb-1"
      const inputClasses =
        value.type === "checkbox" ? "custom-control-input" : "pl-1"
      $("#options-form").append(
        $("<div>", { class: divClasses })
          .append(
            $("<input>", {
              type: value.type,
              class: inputClasses,
              checked: options[key],
              value: value.type === "number" ? options[key] : "",
              max: value.type === "number" ? value.max : "",
              min: value.type === "number" ? value.min : "",
              width: value.type === "number" ? "50px" : "",
            })
              .attr("id", `${key}_input`)
              .change(function() {
                if (value.type == "checkbox") {
                  options[key] = this.checked
                  localStorage.setItem(`${KeyPrefix}_${key}`, this.checked)
                } else {
                  const val = Math.min(
                    Math.max(this.value, value.min),
                    value.max
                  )
                  options[key] = val
                  localStorage.setItem(`${KeyPrefix}_${key}`, val)
                }

                if (key === "showArea") {
                  if (this.checked) {
                    $(".nv-groups").removeClass("area-hidden")
                  } else {
                    $(".nv-groups").addClass("area-hidden")
                  }
                } else if (key === "focusEnable") {
                  chart.options(chartOptions())
                  $(".focusnote").toggle()
                  $(".resetXdomain").toggle()
                } else if (key === "interpolate") {
                  chart.options(chartOptions())
                }

                if (key === "normMergedEx" || key === "calcEff") {
                  updateChart()
                } else if (key === "precision") {
                  data = []
                  svg.datum(data).call(chart)
                  if (chart) {
                    updateChart()
                  }
                } else if (chart) {
                  refreshChart()
                }
              })
          )
          .append(
            $("<label>", {
              for: `${key}_input`,
              class:
                value.type === "checkbox" ? "custom-control-label" : "ml-2",
            }).html(value.msg)
          )
      )
    })

    // initialize chart
    nv.addGraph(
      function() {
        chart = nv.models.lineChart().options(chartOptions())
        chart.lines.duration(0)
        chart.brushExtent(options.startingBrush)
        chart.interactiveLayer.tooltip.valueFormatter(function(d, i) {
          if (d) {
            return `${Math.round(d * 1000) / 10}%`
          }
          return "--"
        })

        // chart sub-models (ie. xAxis, yAxis, etc) when accessed directly
        // return themselves, not the parent chart, so need to chain separately

        // chart.xAxis.axisLabel('Wavelength (nm)');

        chart.yAxis
          .axisLabel("Normalized Ex/Em/Transmission")
          .tickFormat(d3.format("1%"))
          .axisLabelDistance(25)

        svg.datum(data).call(chart)

        chart.focus.dispatch.on("brush", function() {
          updateGlobalGradient()
        })

        chart.dispatch.on("stateChange", function(e) {
          chart.update()
        })

        try {
          const slider = document.getElementById("y-zoom-slider")

          if (slider) {
            noUiSlider.create(slider, {
              start: [1], // 4 handles, starting at...
              behaviour: "tap-drag", // Move handle on tap, bar is draggable
              orientation: "vertical",
              direction: "rtl",
              range: { min: 0.1, max: 1 },
              format: {
                to: function(value) {
                  return Math.round(value * 100) / 100
                },
                from: function(value) {
                  return value
                },
              },
            })

            // update filter settings when user changes slider
            slider.noUiSlider.on("update", function() {
              const m = chart.yDomain()[0]
              chart.yDomain([m, slider.noUiSlider.get()])
              chart.update()
            })
          }
        } catch (error) {
          console.error(error)
        }
      },
      function() {
        resizeYSlider()
        $(window).resize(function() {
          if ($(document).width() < 576) {
            chart.legend.maxKeyLength(14)
          }
          chart.update()
          resizeYSlider()
        })

        $(".resetXdomain").click(function() {
          chart.brushExtent(options.startingBrush)
          refreshChart()
        })
        $(".scale-btns input").change(function() {
          setYscale(this.value)
        })

        $("#undo-scaling").click(function() {
          unscale_all()
        })
        $(window).on("load", function() {
          chart.update()
        })

        setup_from_url(urlParams)

        if (!options.showArea) {
          setTimeout(function() {
            $(".nv-groups").addClass("area-hidden")
          }, 100)
        }
      }
    )

    $("#fluor-select").select2(
      {
        theme: "bootstrap",
        containerCssClass: ":all:",
        width: "auto",
      },
      $("#fluor-select").removeClass("custom-select")
    )

    if (typeof scopespectra !== "undefined") {
      for (let i = 0; i < scopespectra.length; i++) {
        localData[scopespectra[i].slug] = [scopespectra[i]]
      }
    }

    $("#scale-camera").change(function() {
      updateEfficiency()
      refreshChart()
    })

    $(
      ".invswitch, .filter-selector, .fluor-selector, #light-select, #camera-select, scaleToEC_input, scaleToQY_input, #merge-light-exfilter"
    ).change(function() {
      if (options.oneAtaTime && !$(this).hasClass("invswitch")) {
        $(this)
          .closest(".card")
          .find("input")
          .not(this)
          .prop("checked", false)
      }
      updateChart()
    })

    $(".filter-selector, .invswitch").change(function() {
      $("#config-select").val("")
    })

    $("#config-select").change(function() {
      $(".filter-selector:checked").each(function() {
        $(this).prop("checked", false)
      })

      const selected = $(this).find(":selected")
      $.each(selected.data("filters"), function(i, d) {
        $(`#invswitch-filter-${d[0]}`).prop("checked", d[2])
        $(`#filter-${d[1]}-${d[0]}`).prop("checked", true)
      })
      const wave = selected.data("laser")
      const light = selected.data("light")
      const camera = selected.data("camera")
      const comments = selected.data("comments")

      if (wave) {
        $(`#light-select option[value='laser-${wave}']`).prop("selected", true)
      } else if (light) {
        $(`#light-select option[data-id='${light}']`).prop("selected", true)
      } else {
        $("#light-select option[value='']").prop("selected", true)
      }
      if (camera) {
        $(`#camera-select option[data-id='${camera}']`).prop("selected", true)
      }
      if (comments) {
        $("#config_comments").text(comments)
      } else {
        $("#config_comments").text("")
      }
      updateChart()
    })

    $(".switchmodal").click(function(e) {
      e.preventDefault()
      $("#settingsModal").modal("hide")
      $("#embedModal").modal("show")
    })
  })

  function setYscale(type) {
    // type can be log or linear
    chart.lines.duration(200)
    const cd = chart.yDomain()
    const m = cd[0]
    const n = cd[1]
    if (type == "log") {
      options.scale = "log"
      chart.yDomain([0.001, n])
      chart.yScale(d3.scale.log())
      chart.yAxis.tickValues([0.01, 0.033, 0.1, 0.33, 1])
      refreshChart()
    } else {
      options.scale = "linear"
      chart.yDomain([0, n])
      chart.yScale(d3.scale.linear())
      chart.yAxis.tickValues(d3.range(0, 1, 0.2))
      refreshChart()
    }
    chart.lines.duration(0)
  }

  // / Form Events

  let focusedItem
  $("body").on("focus", ".fluor-selector, .select2", function(event) {
    // Store the current value on focus and on change
    focusedItem = $(this)
      .closest(".row")
      .find(".fluor-selector")
      .val()
  })

  // main function when fluor-selector has been changed
  $("body").on("change", ".fluor-selector", function(event) {
    const selector = this
    const slug = $(this).val()
    const row = $(this).closest(".row")
    // if the item is already selected, cancel operation
    if (dataHasSlug(slug) && slug != focusedItem) {
      alert(
        `${localData[slug][0].key
          .replace(" em", "")
          .replace(" 2p", "")
          .replace(" ex", "")} is already selected.`
      )
      if (focusedItem) {
        row
          .find(".fluor-selector")
          .val(focusedItem)
          .change()
      } else {
        row
          .find(".fluor-selector")
          .val(0)
          .change()
      }

      return
    }

    // Remove the previous item fom the list
    removeItem(focusedItem)
    $(selector)
      .siblings(".item-link")
      .remove()

    // Add the new item to the data (unless it's blank)
    // then update the chart
    if (slug != null && slug != 0) {
      addItem(slug).then(function() {
        if (localData[slug][0].url) {
          $(selector)
            .parent()
            .append(
              $("<div>", {
                class: "input-group-append item-link",
                title: "visit item page",
              }).append(
                $("<a>", {
                  class: "item-link",
                  href: localData[slug][0].url,
                  target: "_blank",
                }).append(
                  $("<button>", { class: "btn btn-info", type: "button" }).html(
                    $("#linkSVG").html()
                  )
                )
              )
            )
        }

        // if the slug has a two-2 item...
        if (
          localData[slug]
            .map(function(e) {
              return e.type
            })
            .indexOf(CONST.stype.twop) > -1
        ) {
          // show the checkbox
          row.find("input.2pcheck").prop("disabled", false)
          // get the two photon item in the data
          const item2p = dataItemMatching({
            slug: slug,
            type: CONST.stype.twop,
          })[0] // assume only 1
          // if the 2p button was already checked... keep 2p enabled
          if (row.find("input.2pcheck")[0].checked) {
            item2p.disabled = false
          }
        } else {
          // if it doesn't have 2P data, hide the box...
          row.find("input.2pcheck").prop("disabled", true)
        }
        refreshChart()
      })
    } else {
      refreshChart()
    }

    focusedItem = slug
  })

  $("body").on("click", ".remove-row", function(e) {
    const row = $(this).closest(".row")
    let rowslug = row.find("select").val()
    if (rowslug == "custom_bp" || rowslug == "custom_laser") {
      rowslug = row.attr("id")
    }
    row.remove()
    removeItem(rowslug)

    refreshChart()
    if ($(".fluor-row").length <= 1) {
      $("#toggle_alls").hide()
    }
  })

  // / Form Templates
  function customLaser(wave) {
    return padDataLimits({
      slug: "custom_laser",
      key: `${wave} laser`,
      area: true,
      color: wave_to_color(wave),
      minwave: wave,
      maxwave: wave,
      values: [
        { x: options.precision * Math.round(wave / options.precision), y: 1 },
      ],
    })
  }

  function updateCustomLaser(wave) {
    const F = customLaser(+wave)
    localData[`laser-${wave}`] = [F]
    removeItem("custom_laser")
    pushData(F)
  }

  // / STYLES

  // example:
  // createGradient($('svg')[0],'MyGradient',[
  //   {offset:'5%', 'stop-color':'#f60'},
  //   {offset:'95%','stop-color':'#ff6'}
  // ]);

  const COLORS = {
    380: "#610061",
    381: "#640066",
    382: "#67006a",
    383: "#6a006f",
    384: "#6d0073",
    385: "#6f0077",
    386: "#72007c",
    387: "#740080",
    388: "#760084",
    389: "#780088",
    390: "#79008d",
    391: "#7b0091",
    392: "#7c0095",
    393: "#7c0095",
    394: "#7f009d",
    395: "#8000a1",
    396: "#8100a5",
    397: "#8100a9",
    398: "#8200ad",
    399: "#8200b1",
    400: "#8300b5",
    401: "#8300b9",
    402: "#8300bc",
    403: "#8300c0",
    404: "#8200c4",
    405: "#8200c8",
    406: "#8100cc",
    407: "#8100cf",
    408: "#8000d3",
    409: "#7f00d7",
    410: "#7e00db",
    411: "#7c00de",
    412: "#7b00e2",
    413: "#7900e6",
    414: "#7800e9",
    415: "#7600ed",
    416: "#7400f1",
    417: "#7100f4",
    418: "#6f00f8",
    419: "#6d00fb",
    420: "#6a00ff",
    421: "#6600ff",
    422: "#6100ff",
    423: "#5d00ff",
    424: "#5900ff",
    425: "#5400ff",
    426: "#5000ff",
    427: "#4b00ff",
    428: "#4600ff",
    429: "#4200ff",
    430: "#3d00ff",
    431: "#3800ff",
    432: "#3300ff",
    433: "#2e00ff",
    434: "#2800ff",
    435: "#2300ff",
    436: "#1d00ff",
    437: "#1700ff",
    438: "#1100ff",
    439: "#0a00ff",
    440: "#0000ff",
    441: "#000bff",
    442: "#0013ff",
    443: "#001bff",
    444: "#0022ff",
    445: "#0028ff",
    446: "#002fff",
    447: "#0035ff",
    448: "#003bff",
    449: "#0041ff",
    450: "#0046ff",
    451: "#004cff",
    452: "#0051ff",
    453: "#0057ff",
    454: "#005cff",
    455: "#0061ff",
    456: "#0066ff",
    457: "#006cff",
    458: "#0071ff",
    459: "#0076ff",
    460: "#007bff",
    461: "#007fff",
    462: "#0084ff",
    463: "#0089ff",
    464: "#008eff",
    465: "#0092ff",
    466: "#0097ff",
    467: "#009cff",
    468: "#00a0ff",
    469: "#00a5ff",
    470: "#00a9ff",
    471: "#00aeff",
    472: "#00b2ff",
    473: "#00b7ff",
    474: "#00bbff",
    475: "#00c0ff",
    476: "#00c4ff",
    477: "#00c8ff",
    478: "#00cdff",
    479: "#00d1ff",
    480: "#00d5ff",
    481: "#00daff",
    482: "#00deff",
    483: "#00e2ff",
    484: "#00e6ff",
    485: "#00eaff",
    486: "#00efff",
    487: "#00f3ff",
    488: "#00f7ff",
    489: "#00fbff",
    490: "#00ffff",
    491: "#00fff5",
    492: "#00ffea",
    493: "#00ffe0",
    494: "#00ffd5",
    495: "#00ffcb",
    496: "#00ffc0",
    497: "#00ffb5",
    498: "#00ffa9",
    499: "#00ff9e",
    500: "#00ff92",
    501: "#00ff87",
    502: "#00ff7b",
    503: "#00ff6e",
    504: "#00ff61",
    505: "#00ff54",
    506: "#00ff46",
    507: "#00ff38",
    508: "#00ff28",
    509: "#00ff17",
    510: "#00ff00",
    511: "#09ff00",
    512: "#0fff00",
    513: "#15ff00",
    514: "#1aff00",
    515: "#1fff00",
    516: "#24ff00",
    517: "#28ff00",
    518: "#2dff00",
    519: "#31ff00",
    520: "#36ff00",
    521: "#3aff00",
    522: "#3eff00",
    523: "#42ff00",
    524: "#46ff00",
    525: "#4aff00",
    526: "#4eff00",
    527: "#52ff00",
    528: "#56ff00",
    529: "#5aff00",
    530: "#5eff00",
    531: "#61ff00",
    532: "#65ff00",
    533: "#69ff00",
    534: "#6cff00",
    535: "#70ff00",
    536: "#73ff00",
    537: "#77ff00",
    538: "#7bff00",
    539: "#7eff00",
    540: "#81ff00",
    541: "#85ff00",
    542: "#88ff00",
    543: "#8cff00",
    544: "#8fff00",
    545: "#92ff00",
    546: "#96ff00",
    547: "#99ff00",
    548: "#9cff00",
    549: "#a0ff00",
    550: "#a3ff00",
    551: "#a6ff00",
    552: "#a9ff00",
    553: "#adff00",
    554: "#b0ff00",
    555: "#b3ff00",
    556: "#b6ff00",
    557: "#b9ff00",
    558: "#bdff00",
    559: "#c0ff00",
    560: "#c3ff00",
    561: "#c6ff00",
    562: "#c9ff00",
    563: "#ccff00",
    564: "#cfff00",
    565: "#d2ff00",
    566: "#d5ff00",
    567: "#d8ff00",
    568: "#dbff00",
    569: "#deff00",
    570: "#e1ff00",
    571: "#e4ff00",
    572: "#e7ff00",
    573: "#eaff00",
    574: "#edff00",
    575: "#f0ff00",
    576: "#f3ff00",
    577: "#f6ff00",
    578: "#f9ff00",
    579: "#fcff00",
    580: "#ffff00",
    581: "#fffc00",
    582: "#fff900",
    583: "#fff600",
    584: "#fff200",
    585: "#ffef00",
    586: "#ffec00",
    587: "#ffe900",
    588: "#ffe600",
    589: "#ffe200",
    590: "#ffdf00",
    591: "#ffdc00",
    592: "#ffd900",
    593: "#ffd500",
    594: "#ffd200",
    595: "#ffcf00",
    596: "#ffcb00",
    597: "#ffc800",
    598: "#ffc500",
    599: "#ffc100",
    600: "#ffbe00",
    601: "#ffbb00",
    602: "#ffb700",
    603: "#ffb400",
    604: "#ffb000",
    605: "#ffad00",
    606: "#ffa900",
    607: "#ffa600",
    608: "#ffa200",
    609: "#ff9f00",
    610: "#ff9b00",
    611: "#ff9800",
    612: "#ff9400",
    613: "#ff9100",
    614: "#ff8d00",
    615: "#ff8900",
    616: "#ff8600",
    617: "#ff8200",
    618: "#ff7e00",
    619: "#ff7b00",
    620: "#ff7700",
    621: "#ff7300",
    622: "#ff6f00",
    623: "#ff6b00",
    624: "#ff6700",
    625: "#ff6300",
    626: "#ff5f00",
    627: "#ff5b00",
    628: "#ff5700",
    629: "#ff5300",
    630: "#ff4f00",
    631: "#ff4b00",
    632: "#ff4600",
    633: "#ff4200",
    634: "#ff3e00",
    635: "#ff3900",
    636: "#ff3400",
    637: "#ff3000",
    638: "#ff2b00",
    639: "#ff2600",
    640: "#ff2100",
    641: "#ff1b00",
    642: "#ff1600",
    643: "#ff1000",
    644: "#ff0900",
    645: "#ff0000",
    646: "#ff0000",
    647: "#ff0000",
    648: "#ff0000",
    649: "#ff0000",
    650: "#ff0000",
    651: "#ff0000",
    652: "#ff0000",
    653: "#ff0000",
    654: "#ff0000",
    655: "#ff0000",
    656: "#ff0000",
    657: "#ff0000",
    658: "#ff0000",
    659: "#ff0000",
    660: "#ff0000",
    661: "#ff0000",
    662: "#ff0000",
    663: "#ff0000",
    664: "#ff0000",
    665: "#ff0000",
    666: "#ff0000",
    667: "#ff0000",
    668: "#ff0000",
    669: "#ff0000",
    670: "#ff0000",
    671: "#ff0000",
    672: "#ff0000",
    673: "#ff0000",
    674: "#ff0000",
    675: "#ff0000",
    676: "#ff0000",
    677: "#ff0000",
    678: "#ff0000",
    679: "#ff0000",
    680: "#ff0000",
    681: "#ff0000",
    682: "#ff0000",
    683: "#ff0000",
    684: "#ff0000",
    685: "#ff0000",
    686: "#ff0000",
    687: "#ff0000",
    688: "#ff0000",
    689: "#ff0000",
    690: "#ff0000",
    691: "#ff0000",
    692: "#ff0000",
    693: "#ff0000",
    694: "#ff0000",
    695: "#ff0000",
    696: "#ff0000",
    697: "#ff0000",
    698: "#ff0000",
    699: "#ff0000",
    700: "#ff0000",
    701: "#fd0000",
    702: "#fb0000",
    703: "#fa0000",
    704: "#f80000",
    705: "#f60000",
    706: "#f40000",
    707: "#f20000",
    708: "#f10000",
    709: "#ef0000",
    710: "#ed0000",
    711: "#eb0000",
    712: "#e90000",
    713: "#e80000",
    714: "#e60000",
    715: "#e40000",
    716: "#e20000",
    717: "#e00000",
    718: "#de0000",
    719: "#dc0000",
    720: "#db0000",
    721: "#d90000",
    722: "#d70000",
    723: "#d50000",
    724: "#d30000",
    725: "#d10000",
    726: "#cf0000",
    727: "#ce0000",
    728: "#cc0000",
    729: "#ca0000",
    730: "#c80000",
    731: "#c60000",
    732: "#c40000",
    733: "#c20000",
    734: "#c00000",
    735: "#be0000",
    736: "#bc0000",
    737: "#ba0000",
    738: "#b90000",
    739: "#b70000",
    740: "#b50000",
    741: "#b30000",
    742: "#b10000",
    743: "#af0000",
    744: "#ad0000",
    745: "#ab0000",
    746: "#a90000",
    747: "#a70000",
    748: "#a50000",
    749: "#a30000",
    750: "#a10000",
    751: "#9f0000",
    752: "#9d0000",
    753: "#9b0000",
    754: "#990000",
    755: "#970000",
    756: "#950000",
    757: "#930000",
    758: "#910000",
    759: "#8f0000",
    760: "#8d0000",
    761: "#8a0000",
    762: "#880000",
    763: "#860000",
    764: "#840000",
    765: "#820000",
    766: "#800000",
    767: "#7e0000",
    768: "#7c0000",
    769: "#7a0000",
    770: "#770000",
    771: "#750000",
    772: "#730000",
    773: "#710000",
    774: "#6f0000",
    775: "#6d0000",
    776: "#6a0000",
    777: "#680000",
    778: "#660000",
    779: "#640000",
    780: "#610000",
  }

  function wave_to_color(wave) {
    wave = Math.round(wave)
    if (wave < 380) {
      return "#000061"
    }
    if (wave > 780) {
      return "#610000"
    }
    return COLORS[wave]
  }

  const chartsvg = $("#spectra svg")[0]
  const chartsvgNS = chartsvg.namespaceURI
  const grad = document.createElementNS(chartsvgNS, "linearGradient")
  grad.setAttribute("id", "wavecolor_gradient")
  grad.setAttribute("class", "svggradient")
  const defs =
    chartsvg.querySelector("defs") ||
    chartsvg.insertBefore(
      document.createElementNS(chartsvgNS, "defs"),
      svg.firstChild
    )
  defs.appendChild(grad)

  function updateGlobalGradient() {
    $(grad).empty()
    const cmin = chart.xAxis.domain()[0]
    const cmax = chart.xAxis.domain()[1]
    const range = cmax - cmin
    let stop
    for (let w = cmin; w < cmax; w += 50) {
      stop = document.createElementNS(chartsvgNS, "stop")
      stop.setAttribute("offset", `${Math.round((100 * (w - cmin)) / range)}%`)
      stop.setAttribute("stop-color", wave_to_color(w))
      grad.appendChild(stop)
    }
    stop = document.createElementNS(chartsvgNS, "stop")
    stop.setAttribute("offset", "100%")
    stop.setAttribute("stop-color", wave_to_color(cmax))
    grad.appendChild(stop)
  }

  // // EFFICIENCY CALCULATIONS

  function spectral_product(ar1, ar2) {
    // calculate product of two spectra.values
    const output = []
    const step = ar1[1].x - ar1[0].x
    const left = Math.max(ar1[0].x, ar2[0].x)
    const right = Math.min(ar1[ar1.length - 1].x, ar2[ar2.length - 1].x)

    const a1 = ar1.slice(
      ar1.findIndex(i => i.x === left),
      ar1.findIndex(i => i.x === right)
    )
    const a2 = ar2.slice(
      ar2.findIndex(i => i.x === left),
      ar2.findIndex(i => i.x === right)
    )
    for (let i = 0; i < a1.length; i++) {
      output.push({ x: a1[i].x, y: a1[i].y * a2[i].y })
    }

    // var offsetA1 = (left - ar1[0].x) / step; // these assume monotonic increase w/ step = 1
    // var offsetA2 = (left - ar2[0].x) / step; // these assume monotonic increase w/ step = 1
    // for (var i = 0; i < right - left; i++) {
    //     console.log(ar1[offsetA1 + i].x, ar2[offsetA2 + i])
    //     if (ar1[offsetA1 + i] && ar2[offsetA2 + i]){
    //         output.push({ x: ar1[offsetA1 + i].x, y: ar1[offsetA1 + i].y * ar2[offsetA2 + i].y });
    //     }
    // }
    return output
  }

  function trapz(arr, min, max) {
    // approximate area under curve as series of trapezoids
    min = min || 300
    max = max || 1000
    let sum = 0
    const step = arr[1].x - arr[0].x
    for (let i = 0; i < arr.length - 1; i++) {
      if (arr[i].x > min) {
        const d = (step * (arr[i].y + arr[i + 1].y)) / 2
        sum += d
      }
      if (arr[i].x > max) {
        break
      }
    }
    return sum
  }

  function invertData(values) {
    const out = []
    for (let i = 0; i < values.length; i++) {
      out.push({ x: values[i].x, y: 1 - values[i].y })
    }
    return out
  }

  function normalizeSpectrum(specvals, maxY) {
    maxY =
      maxY ||
      Math.max.apply(
        Math,
        specvals.map(function(a) {
          return a.y
        })
      )
    for (let i = 0; i < specvals.length; i++) {
      specvals[i].y /= maxY
    }
    return [specvals, maxY]
  }

  function combineSpectra(pathlist) {
    return pathlist.reduce(function(acc, cur) {
      if (acc) {
        return spectral_product(acc, cur.values)
      }
      return cur.values
    }, null)
  }

  // disgusting awful code design
  function calcExEmPaths() {
    const empath = []
    const expath = []
    const bspath = $(".bs-filter:checked")
      .map(function() {
        return data.find(d => d.slug === this.value)
      })
      .get()

    for (let n = 0; n < bspath.length; n++) {
      const invVals = { ...bspath[n] }
      invVals.values = invertData(bspath[n].values)
      empath.push(bspath[n])
      expath.push(invVals)
    }

    empath.push.apply(
      empath,
      $(".em-filter:checked")
        .map(function() {
          return data.find(d => d.slug === this.value)
        })
        .get()
    )
    if ($("#scale-camera").prop("checked")) {
      if ($("#camera-select :selected").val()) {
        empath.push(
          data.find(d => d.slug === $("#camera-select :selected").val())
        )
      }
    }

    if (
      $("#merge-light-exfilter").prop("checked") &&
      $("#light-select :selected").val() !== ""
    ) {
      const d = data.find(d => d.slug.indexOf("merged") >= 0)
      if (d) {
        expath.push(d)
      }
    } else {
      expath.push.apply(
        expath,
        $(".ex-filter:checked")
          .map(function() {
            return data.find(d => d.slug === this.value)
          })
          .get()
      )
      if ($("#light-select :selected").val()) {
        expath.push(
          data.find(d => d.slug === $("#light-select :selected").val())
        )
      }
    }
    return {
      exvalues: combineSpectra(expath),
      emvalues: combineSpectra(empath),
    }
  }

  function efficiency(filtervalues, fluor) {
    const dataitem = dataItemMatching(fluor)[0]
    if (!dataitem) return
    const fluorspectrum = dataitem.values
    const filter_dye_spectrum = spectral_product(filtervalues, fluorspectrum)
    const filter_dye_area = trapz(filter_dye_spectrum)
    if (fluor.type == "em") {
      var efficiency = filter_dye_area / trapz(fluorspectrum)
    } else if (fluor.type == "ex") {
      if (options.exEffBroadband) {
        var efficiency = filter_dye_area / trapz(fluorspectrum)
      } else {
        var efficiency = filter_dye_area / trapz(filtervalues)
      }
    }

    return {
      key: `${dataitem.key} eff`,
      values: filter_dye_spectrum,
      efficiency: efficiency,
      area: true,
      color:
        fluor.type == "em"
          ? "url(#diagonal-stripe-r)"
          : "url(#diagonal-stripe-l)",
      classed: `${fluor.type}-efficiency subtype-eff`,
      type: "eff",
      scalar: dataitem.scalar,
    }
  }

  d3.selection.prototype.moveToFront = function() {
    return this.each(function() {
      this.parentNode.appendChild(this)
    })
  }

  d3.selection.prototype.moveToBack = function() {
    return this.each(function() {
      const { firstChild } = this.parentNode
      if (firstChild) {
        this.parentNode.insertBefore(this, firstChild)
      }
    })
  }
  Array.prototype.clean = function(deleteValue) {
    for (let i = 0; i < this.length; i++) {
      if (this[i] == deleteValue) {
        this.splice(i, 1)
        i--
      }
    }
    return this
  }

  function autoSizeText() {
    let el
    let elements
    let _i
    let _len
    let _results
    elements = $(".filter-label")
    if (elements.length < 0) {
      return
    }
    _results = []
    for (_i = 0, _len = elements.length; _i < _len; _i++) {
      el = elements[_i]
      _results.push(
        (function(el) {
          let resizeText
          let _results1
          resizeText = function() {
            let elNewFontSize
            elNewFontSize = `${parseInt(
              $(el)
                .css("font-size")
                .slice(0, -2)
            ) - 0.5}px`
            return $(el).css("font-size", elNewFontSize)
          }
          _results1 = []
          let _c = 0
          while ((el.scrollHeight > el.offsetHeight) & (_c < 20)) {
            _results1.push(resizeText())
            _c += 1
          }
          return _results1
        })(el)
      )
    }
    return _results
  }

  function selectedSlugs() {
    const selected = $(".filter-selector:checked, .data-selector :selected")
    let list
    if (
      $("#merge-light-exfilter").prop("checked") &&
      $("#light-select :selected").val() !== ""
    ) {
      list = selected
        .not("#light-select :selected, .ex-filter:checked")
        .map(function() {
          return $(this).val()
        })
      let _concat = "merged/"
      _concat += $("#light-select :selected, .ex-filter:checked")
        .map(function(i, o) {
          return o.value
        })
        .get()
        .reduce(function(a, b) {
          return `${a}/${b}`
        })
      if (options.normMergedEx) {
        _concat += "/normed"
      }
      list.push(_concat)
    } else {
      list = selected.map(function() {
        return $(this).val()
      })
    }
    return list.get().clean("")
  }

  function updateData() {
    const selected = selectedSlugs()
    // remove old data that isnt still in the array
    for (let i = data.length - 1; i >= 0; i--) {
      if ($.inArray(data[i].slug, selected) == -1) {
        data.splice(i, 1)
      }
    }
    // add new items
    const deferreds = []
    for (let n = 0; n < selected.length; n++) {
      if (selected[n] && !dataHasSlug(selected[n])) {
        if (selected[n].startsWith("laser")) {
          updateCustomLaser(selected[n].replace("laser-", ""))
        } else if (selected[n].startsWith("merged")) {
          deferreds.push(addMerged(selected[n]))
        } else {
          deferreds.push(addItem(selected[n]))
        }
      }
    }
    return deferreds
  }

  function invertInverted() {
    const inverted = $(".invswitch:checked")
      .map(function() {
        return $(this).val()
      })
      .get()
    for (let n = 0; n < data.length; n++) {
      if (data[n].inverted & ($.inArray(data[n].slug, inverted) == -1)) {
        data[n].values = invertData(data[n].values)
        data[n].inverted = false
      } else if (
        (typeof data[n].inverted === "undefined" || !data[n].inverted) &
        ($.inArray(data[n].slug, inverted) > -1)
      ) {
        data[n].values = invertData(data[n].values)
        data[n].inverted = true
      }
    }
  }

  function updateChart() {
    const deferreds = updateData()
    if (deferreds) {
      $.when.apply($, deferreds).then(function() {
        invertInverted()
        updateEfficiency()
        refreshChart()
        d3.selectAll(".subtype-pd").moveToBack()
        d3.selectAll(".subtype-qe").moveToBack()
        d3.selectAll(".subtype-eff").moveToFront()
      })
    } else {
      updateEfficiency()
    }
  }

  function updateEfficiency() {
    if (!options.calcEff) {
      $("#efficiency-text").attr("placeholder", "disabled")
      return
    }

    const a = calcExEmPaths()
    const L = $("#light-select :selected")
    if (L.data("type") === "laser") {
      const Lwave = +L.val().replace("laser-", "")
      try {
        var emeff = a.emvalues.filter(function(d) {
          return d.x === Lwave || d.x === Lwave + 1
        })[0].y
      } catch (e) {
        var emeff = null
      }
      var cameff = 1
      try {
        if ($("#scale-camera").prop("checked")) {
          var cameff = dataItemMatching({ category: "c" })[0].values.filter(
            function(d) {
              return d.x === Lwave || d.x === Lwave + 1
            }
          )[0].y
        }
      } catch (e) {}
      if (emeff) {
        const blocking = Math.round(-Math.log10(emeff / cameff) * 10) / 10
        $("#laserBlocking").text(`OD ${blocking}`)
      }
    } else {
      $("#laserBlocking").text("")
    }

    if (!$("#fluor-select :selected").val()) {
      return
    }
    removeSubtypes("eff")
    const fluorslug = $("#fluor-select :selected").val()
    if (a.emvalues && a.emvalues.length) {
      var emeff = efficiency(a.emvalues, { slug: fluorslug, type: "em" })
    }
    if (a.exvalues && a.exvalues.length) {
      var exeff = efficiency(a.exvalues, { slug: fluorslug, type: "ex" })
    }
    let info = ""
    if (exeff) {
      // if (options.normMergedEx)
      // exeff.values = normalizeSpectrum(exeff.values, options.normMergedScalar)[0]
      pushData(exeff)
      if (exeff.efficiency)
        info += `Ex: ${String(Math.round(exeff.efficiency * 1000) / 10)}%`
    }
    if (emeff) {
      pushData(emeff)
      if (emeff.efficiency) {
        if (info) info += ", "
        info += `Em: ${String(Math.round(emeff.efficiency * 1000) / 10)}%`
      }
    }

    if ((emeff != undefined) & (exeff != undefined)) {
      const b =
        exeff.efficiency * exeff.scalar * emeff.efficiency * emeff.scalar
      if (b) info = `${String(Math.round(b / 10) / 100)} (${info})`
    }

    if (info != "") {
      $("#efficiency-text").attr("placeholder", info)
    } else {
      $("#efficiency-text").attr("placeholder", "Select a fluor and config...")
    }
  }

  $(window).on("load", function() {
    setTimeout(function() {
      autoSizeText()
      chart && chart.update()
      if ($(document).width() < 576) {
        chart.legend.maxKeyLength(15)
      }

      $(".filter-label").each(function(i) {
        this.innerHTML = this.innerHTML.replace(/^\s*\w+/, function(x) {
          const f = x.trim().toLowerCase()
          if (
            f === "chroma" ||
            f === "semrock" ||
            f === "zeiss" ||
            f === "omega"
          ) {
            return `<span class="first-word">${x}</span>`
          }
          return x
        })
      })
    }, 150)
  })

  if ($(".microscope-wrapper").length) {
    const topofDiv = $(".microscope-wrapper").offset().top
    $(window).scroll(function() {
      if (options.stickySpectra) {
        if ($(window).scrollTop() > topofDiv) {
          $(".microscope-wrapper").addClass("shadowed")
        } else {
          $(".microscope-wrapper").removeClass("shadowed")
        }
      }
    })
  }

  $("#stickyPin").click(function() {
    $(".pin-wrapper").toggleClass("rotate-90")
    if ($(".pin-wrapper").hasClass("rotate-90")) {
      options.stickySpectra = false
      $(".microscope-wrapper").removeClass("sticky")
      $(".microscope-wrapper").removeClass("shadowed")
    } else {
      options.stickySpectra = true
      $(".microscope-wrapper").addClass("sticky")
    }
    localStorage.setItem(`${KeyPrefix}_stickySpectra`, options.stickySpectra)
  })

  function build_current_uri() {
    const params = []
    if ($("#config-select").val()) {
      params.push(`c=${encodeURIComponent($("#config-select").val())}`)
    } else {
      const filters = $(".filter-selector:checked")
        .map(function(i, d) {
          let v = $(d).prop("id")
          if (
            $(this)
              .closest("li")
              .find(".invswitch")
              .prop("checked")
          ) {
            v += "!"
          }
          return encodeURIComponent(v)
        })
        .toArray()
      if (filters.length) {
        params.push(`f=${filters.join(",")}`)
      }
    }
    if ($(".fluor-selector").val()) {
      params.push(`p=${encodeURIComponent($(".fluor-selector").val())}`)
    }
    if ($("#light-select").val()) {
      params.push(`l=${encodeURIComponent($("#light-select").val())}`)
    }
    if ($("#camera-select").val()) {
      params.push(`d=${encodeURIComponent($("#camera-select").val())}`)
    }
    let outurl = `${window.location.protocol}//${window.location.host}${window.location.pathname}`
    if (params.length) {
      outurl += `?${params.join("&")}`
    }
    return outurl
  }

  function copyToClipboard(elem) {
    // create hidden text element, if it doesn't already exist
    const targetId = "_hiddenCopyText_"
    const isInput = elem.tagName === "INPUT" || elem.tagName === "TEXTAREA"
    let origSelectionStart
    let origSelectionEnd
    if (isInput) {
      // can just use the original source element for the selection and copy
      target = elem
      origSelectionStart = elem.selectionStart
      origSelectionEnd = elem.selectionEnd
    } else {
      // must use a temporary form element for the selection and copy
      target = document.getElementById(targetId)
      if (!target) {
        var target = document.createElement("textarea")
        target.style.position = "absolute"
        target.style.left = "-9999px"
        target.style.top = "0"
        target.id = targetId
        document.body.appendChild(target)
      }
      target.textContent = elem.textContent
    }
    // select the content
    const currentFocus = document.activeElement
    target.focus()
    target.setSelectionRange(0, target.value.length)

    // copy the selection
    let succeed
    try {
      succeed = document.execCommand("copy")
    } catch (e) {
      succeed = false
    }
    // restore original focus
    if (currentFocus && typeof currentFocus.focus === "function") {
      currentFocus.focus()
    }

    if (isInput) {
      // restore prior selection
      elem.setSelectionRange(origSelectionStart, origSelectionEnd)
    } else {
      // clear temporary content
      target.textContent = ""
    }
    return succeed
  }

  $("#copyButton").click(function() {
    $("#shareLink").prop("disabled", false)
    copyToClipboard(document.getElementById("shareLink"))
    $("#shareLink").prop("disabled", true)
  })

  $("#share-button").click(function() {
    const link = build_current_uri()
    $("#shareLink").val(link)
    $("#mailLink").attr(
      "href",
      $("#mailLink").data("mailinfo") + encodeURIComponent(link)
    )
    $("#shareModal").modal("show")
  })
})()

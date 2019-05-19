/* eslint-disable no-plusplus */
import d3 from "d3";
import noUiSlider from "nouislider";
import $ from "jquery";
import nv from "nvd3";
import {
  options,
  userOptions,
  CONST
} from "./constants";
import {
  formSelection,
  fluorRow,
  excRow,
  filterRow
} from "./spectra_form";
import COLORS from "./colors";

export default function initSpectra(selection, _spectra_options) {
  let spectra_options = _spectra_options
  let chart;
  const data = [];
  const localData = {};
  const svg = d3.select(selection);
  let overlaps = [];
  let disabledData = [];

  function padDataLimits(d) {
    for (let i = d.minwave - 1; i >= options.minwave; i--) {
      d.values.unshift({
        x: i,
        y: 0
      });
    }
    for (let n = d.maxwave + 1; n <= Math.max(options.maxwave, 1000); n++) {
      d.values.push({
        x: n,
        y: 0
      });
    }
    return d;
  }

  function getData(slug) {
    const dfd = $.Deferred();
    // download if not already downloaded
    if (!(slug in localData)) {
      $.getJSON(slug)
        .done(d => {
          for (let n = 0; n < d.spectra.length; n++) {
            d.spectra[n] = padDataLimits(d.spectra[n]);
            d.spectra[n].exNormed = 1;
            if (d.spectra[n].type === CONST.stype.twop && options.hide2p) {
              d.spectra[n].disabled = true;
            } else if (
              d.spectra[n].type === CONST.stype.twop &&
              !options.hide2p
            ) {
              d.spectra[n].disabled = false;
            }
            if (
              d.spectra[n].category === CONST.category.light ||
              d.spectra[n].category === CONST.category.camera
            ) {
              d.spectra[n].gradient = true;
            }
          }
          localData[slug] = d.spectra;
          dfd.resolve(localData[slug]);
        })
        .fail(d => {
          dfd.reject(d.status);
        });
    } else {
      // otherwise pull from local dict
      dfd.resolve(localData[slug]);
    }
    return dfd.promise();
  }

  function dataHasKey(key) {
    return $.grep(data, obj => obj.key === key).length > 0;
  }

  function dataHasSlug(slug) {
    return $.grep(data, obj => obj.slug === slug).length > 0;
  }

  function slugOptions(slug) {
    for (let i = 0; i < spectra_options.length; i++) {
      if (spectra_options[i].slug == slug) {
        return spectra_options[i];
      }
    }
    return false;
  }

  function dataItemMatching(filter, d) {
    d = d || data;
    return d.filter(item => {
      for (const key in filter) {
        if (typeof item[key] === "undefined" || item[key] != filter[key])
          return false;
      }
      return true;
    });
  }

  function addItem(slug, subtype_) {
    // add spectral data to array
    const subtype = subtype_ || false;

    return getData(slug).then(d => {
      for (let i = 0; i < d.length; i++) {
        if ((!subtype || d[i].type === subtype) && !dataHasKey(d[i].key)) {
          data.push(JSON.parse(JSON.stringify(d[i]))); // make a copy of the object
        }
      }
    });
  }

  function removeItem(slug, subtype_) {
    // optional subtype argument can be used to remove a specific type of spectrum
    // rather than the whole slug family
    const subtype = subtype_ || false;

    const ix = $.map(data, (obj, index) => {
      if (obj.slug === slug) {
        if (!subtype || obj.type === subtype) {
          return index;
        }
      }
    });
    for (let i = ix.length - 1; i >= 0; i--) data.splice(ix[i], 1);
  }

  // /// chart setup

  function unscaleAll() {
    options.exNormWave = undefined;
    options.scaleToQY = false;
    options.scaleToEC = false;
    $("#exnormRadioOFF").prop("checked", true);
    $("#scaleToQY_input").prop("checked", false);
    $("#scaleToEC_input").prop("checked", false);
    refreshChart();
  }

  function scaleDataUp(filter) {
    // scale data "up" according to the data.scalar value
    // filter can be .e.g. {type: 'ex'} to scale by ExtCoeff
    const maxScalar =
      Math.max.apply(
        null,
        data.map(function (e) {
          return e.scalar || 0;
        })
      ) || 1;
    for (let n = 0; n < data.length; n++) {
      // only scale certain data by filter
      let skip = false;
      if (data[n].scaled || typeof data[n].scalar === "undefined") {
        skip = true;
      }
      for (const key in filter) {
        if (
          typeof data[n][key] === "undefined" ||
          data[n][key] != filter[key]
        ) {
          skip = true;
          break;
        }
      }

      if (!skip) {
        let SCALE = data[n].scalar || 0.001;
        if (data[n].type == "ex") {
          SCALE /= maxScalar;
        }
        // do the scaling
        for (let i = 0; i < data[n].values.length; i++) {
          data[n].values[i].y *= SCALE;
        }
        data[n].scaled = SCALE;
      }
    }
  }

  function scaleDataDown(filter) {
    for (let n = 0; n < data.length; n++) {
      // only scale certain data by filter
      let skip = false;
      for (const key in filter) {
        if (
          typeof data[n][key] === "undefined" ||
          data[n][key] != filter[key] ||
          !data[n].scaled
        ) {
          skip = true;
          break;
        }
      }
      if (!skip) {
        for (let i = 0; i < data[n].values.length; i++) {
          data[n].values[i].y /= data[n].scaled;
        }
        data[n].scaled = false;
      }
    }
  }

  function excitationNormalization() {
    if (options.exNormWave) {
      for (let n = 0; n < data.length; n++) {
        if (data[n].type === "em") {
          const exvals = dataItemMatching({
              slug: data[n].slug,
              type: "ex"
            })[0]
            .values;
          let targetScalar = dataItemMatching({
              x: options.exNormWave
            },
            exvals
          )[0].y;
          targetScalar = Math.max(targetScalar, 0.0001); // so we don't clobber the values with zero
          const currentScalar = data[n].exNormed;
          if (currentScalar !== targetScalar) {
            const S = targetScalar / currentScalar;
            data[n].values = data[n].values.map(item => {
              return {
                x: item.x,
                y: item.y * S
              };
            });
            data[n].exNormed = targetScalar;
          }
        }
      }
    } else {
      for (let n = 0; n < data.length; n++) {
        if (data[n].type === "em" && data[n].exNormed !== 1) {
          data[n].values = data[n].values.map(item => {
            return {
              x: item.x,
              y: item.y / data[n].exNormed
            };
          });
          data[n].exNormed = 1;
        }
      }
    }
  }

  function resizeYSlider() {
    try {
      const rect = $("#spectrasvg")
        .find(".nv-focus")
        .first()
        .find(".nv-background>rect");
      const h = +rect.attr("height");
      // var t = Math.abs(rect.offset().top - $('#spectrasvg').offset().top);
      $("#y-zoom-slider").height(h * 0.88);
    } catch (e) {}
  }

  function scaleDataToOptions() {
    if (options.scaleToEC) {
      scaleDataUp({
        type: "ex"
      });
    } else {
      scaleDataDown({
        type: "ex"
      });
    }

    if (options.scaleToQY) {
      scaleDataUp({
        type: "em"
      });
    } else {
      scaleDataDown({
        type: "em"
      });
    }
  }

  function refreshChart() {
    chart.lines.duration(300);
    if (options.autoscaleBrush) {
      let smin = 10000;
      let smax = 0;
      for (let i = 0; i < data.length; i++) {
        if (!data[i].disabled) {
          smin = Math.min(data[i].minwave, smin);
          smax = Math.max(data[i].maxwave, smax);
        }
      }
      chart.brushExtent([smin, smax]);
    }

    scaleDataToOptions();
    excitationNormalization();
    calculateEfficiency();
    chart.update();
    updateGlobalGradient();
    chart.lines.duration(0);
    if (options.scaleToQY || options.scaleToEC || options.exNormWave) {
      $("#y-zoom-slider").show();
      resizeYSlider();
    } else {
      $("#y-zoom-slider").hide();
    }
  }

  function addFormItem(category, stype, open_, value_) {
    const open = open_ || false;
    const value = value_ || undefined;
    const filter = {
      category
    };
    if (stype) {
      filter.subtype = stype;
    }
    const selWidget = formSelection(filter, spectra_options);

    if (category === CONST.category.protein) {
      $(fluorRow(selWidget)).appendTo($("#protein-table"));
    }
    if (category === CONST.category.dye) {
      $(fluorRow(selWidget)).appendTo($("#dye-table"));
    }
    if (
      category === CONST.category.dye ||
      category === CONST.category.protein
    ) {
      if ($(".fluor-row").length > 1) {
        $("#toggle_alls").show();
      }
    } else if (category == CONST.category.light) {
      $(excRow(selWidget, "light-row", refreshChart)).appendTo(
        $("#light-table")
      );
    } else if (
      category === CONST.category.filter &&
      (stype === CONST.stype.bpx || stype === CONST.stype.sp)
    ) {
      $(filterRow(selWidget, stype)).appendTo($("#exfilter-table"));
    } else if (
      category === CONST.category.filter &&
      (stype === CONST.stype.bpm || stype === CONST.stype.lp)
    ) {
      $(filterRow(selWidget, stype)).appendTo($("#emfilter-table"));
    } else if (category === CONST.category.camera) {
      $(filterRow(selWidget, stype)).appendTo($("#camqe-table"));
    }

    const a = selWidget.select2({
      theme: "bootstrap"
    });
    if (value) {
      a.val(value).change();
    } else if (open) {
      focusedItem = $(this)
        .closest(".row")
        .find(".data-selector")
        .val();
      a.select2("open");
    }
  }

  // // Scaling functions

  // // ON WINDOW LOAD
  function getUrlParams(prop) {
    const params = {};
    const search = decodeURIComponent(
      window.location.href.slice(window.location.href.indexOf("?") + 1)
    );
    const definitions = search.split("&");

    definitions.forEach(function (val, key) {
      const parts = val.split("=", 2);
      params[parts[0]] = parts[1];
    });

    return prop && prop in params ? params[prop] : params;
  }

  $(function pageLoad() {
    $("#y-zoom-slider").hide();
    $('[data-toggle="popover"]').popover();

    $.each(userOptions, function setUpUserOptions(key, value) {
      $("#options-form").append(
        $("<div>", {
          class: "custom-control custom-checkbox mb-1 pb-1"
        })
        .append(
          $("<input>", {
            type: value.type,
            class: "custom-control-input",
            checked: options[key]
          })
          .attr("id", `${key}_input`)
          .change(e => {
            if (value.type === "checkbox") {
              options[key] = e.target.checked;
            } else {
              options[key] = e.target.value;
            }

            if (key === "showArea") {
              if (e.target.checked) {
                $(".nv-groups").removeClass("area-hidden");
              } else {
                $(".nv-groups").addClass("area-hidden");
              }
            }
            refreshChart();
          })
        )
        .append(
          $("<label>", {
            for: `${key}_input`,
            class: "custom-control-label"
          }).text(value.msg)
        )
      );
    });

    // initialize chart
    nv.addGraph(() => {
      chart = nv.models.lineChart().options({
        focusEnable: true,
        focusShowAxisX: false,
        noData: "Add spectra below ...",
        showLegend: true,
        showXAxis: true,
        showYAxis: true,
        duration: 300,
        useInteractiveGuideline: !window.mobilecheck(),
        clipEdge: true,
        margin: {
          left: 40,
          bottom: 15
        },
        yDomain: [0, 1]
        // forceY: [0,1.04],
        // forceX: [350, 750],
      });
      chart.lines.duration(0);
      chart.brushExtent(options.startingBrush);
      chart.interactiveLayer.tooltip.valueFormatter(d => {
        const out = d ? `${Math.round(d * 1000) / 10}%` : "--";
        return out;
      });

      // chart sub-models (ie. xAxis, yAxis, etc) when accessed directly
      // return themselves, not the parent chart, so need to chain separately

      // chart.xAxis.axisLabel('Wavelength (nm)');

      chart.yAxis
        .axisLabel("Normalized Ex/Em/Transmission")
        .tickFormat(d3.format("1%"))
        .axisLabelDistance(25);

      svg.datum(data).call(chart);

      chart.focus.dispatch.on("brush", () => {
        updateGlobalGradient();
      });

      chart.dispatch.on("stateChange", () => {
        chart.update();
      });

      const slider = document.getElementById("y-zoom-slider");
      noUiSlider.create(slider, {
        start: [1], // 4 handles, starting at...
        behaviour: "tap-drag", // Move handle on tap, bar is draggable
        orientation: "vertical",
        direction: "rtl",
        range: {
          min: 0.1,
          max: 1
        },
        format: {
          to(value) {
            return Math.round(value * 100) / 100;
          },
          from(value) {
            return value;
          }
        }
      });

      // update filter settings when user changes slider
      slider.noUiSlider.on("update", () => {
        const m = chart.yDomain()[0];
        chart.yDomain([m, slider.noUiSlider.get()]);
        chart.update();
      });
    });

    resizeYSlider();
    $(window).resize(() => {
      chart.update();
      resizeYSlider();
    });

    const urlParams = getUrlParams();
    if ("s" in urlParams) {
      const arr = urlParams.s.toLowerCase().split(",");
      for (let i = 0; i < arr.length; i++) {
        const opts = slugOptions(arr[i]);
        if (opts) {
          try {
            // STILL BUGGY
            addFormItem(opts.category, opts.subtype, false, opts.slug);
          } catch (e) {
            console.log(e);
          }
        }
      }
    } else {
      addFormItem("p", null, false, "egfp_default");
      addFormItem("d");
    }
    addFormItem("l");
    addFormItem("f", "bx");
    addFormItem("f", "bm");
    addFormItem("c");

    $(".scale-btns input").change(e => setYscale(e.target.value));

    $("body").on("change", ".singlecheck", function (e) {
      // when the ex/em/2p checkbox is checked,
      // enable/disable the corresponding item in the data
      const slug = $(this)
        .closest(".row")
        .find("select")
        .val();
      const type = $(this).val();
      for (let i = 0; i < data.length; i++) {
        if (data[i].slug == slug && data[i].type == type) {
          data[i].disabled = !this.checked;
        }
        if (type == CONST.stype.twop && localData[slug]) {
          const pos = localData[slug]
            .map(function (e) {
              return e.type;
            })
            .indexOf("2p");
          if (localData[slug][pos]) {
            localData[slug][pos].disabled = !this.checked;
          }
        }
      }
      refreshChart();
    });

    $(".toggleall").change(function (e) {
      $(`.${$(this).val()}check`)
        .prop("checked", $(this).is(":checked"))
        .change();
    });

    $(".singlecheck").change(function () {
      const c = $(`.singlecheck.${$(this).val()}check`).map(function () {
        return $(this).is(":checked");
      });
      if (
        c.toArray().every(function (a) {
          return a;
        })
      ) {
        $(`#toggle_all_${$(this).val()}`).prop("indeterminate", false);
        $(`#toggle_all_${$(this).val()}`).prop("checked", true);
      } else if (
        c.toArray().every(function (a) {
          return !a;
        })
      ) {
        $(`#toggle_all_${$(this).val()}`).prop("indeterminate", false);
        $(`#toggle_all_${$(this).val()}`).prop("checked", false);
      } else {
        $(`#toggle_all_${$(this).val()}`).prop("indeterminate", true);
      }
    });

    $("body").on("click", ".nv-focusWrap", function () {
      // if the user moves the focus, don't autoscale on them
      options.autoscaleBrush = false;
      $('#options_form input[data-opt="autoscaleBrush"]').prop(
        "checked",
        false
      );
    });

    // eyeSVG = $("#eyeSVG").html();
    // linkSVG = $("#linkSVG").html();

    $("#undo-scaling").click(function () {
      unscaleAll();
    });

    $(".resetXdomain").click(function () {
      chart.brushExtent(options.startingBrush);
      refreshChart();
    });
  });

  function setYscale(type) {
    // type can be log or linear
    chart.lines.duration(300);
    const cd = chart.yDomain();
    const m = cd[0];
    const n = cd[1];
    if (type == "log") {
      options.scale = "log";
      chart.yDomain([0.001, n]);
      chart.yScale(d3.scale.log());
      chart.yAxis.tickValues([0.01, 0.033, 0.1, 0.33, 1]);
      refreshChart();
    } else {
      options.scale = "linear";
      chart.yDomain([0, n]);
      chart.yScale(d3.scale.linear());
      chart.yAxis.tickValues(d3.range(0, 1, 0.2));
      refreshChart();
    }
    chart.lines.duration(0);
  }

  function waveToColor(wave_) {
    const wave = Math.round(wave_);
    if (wave < 380) {
      return "#000061";
    }
    if (wave > 780) {
      return "#610000";
    }
    return COLORS[wave];
  }

  // / Form Templates
  function customLaser(row) {
    const id = row.attr("id");
    const wave = +row.find('input[name="custom_laser_wave"]').val();
    return padDataLimits({
      slug: id,
      key: `${wave} ex`,
      area: true,
      minwave: wave,
      maxwave: wave,
      values: [{
        x: wave,
        y: 1
      }]
    });
  }

  function customFilter(row) {
    const id = row.attr("id");
    const center = row.find('input[name="custom_em_center"]').val();
    const width = row.find('input[name="custom_em_bandwidth"]').val();
    const trans = row.find('input[name="custom_em_trans"]').val();
    const minX = parseFloat(center) - width / 2;
    const maxX = parseFloat(center) + width / 2;
    const vals = [];
    for (let n = Math.min(minX, 300); n < Math.max(maxX, 1000); n++) {
      if (n >= minX && n <= maxX) {
        vals.push({
          x: n,
          y: +trans
        });
      } else {
        vals.push({
          x: n,
          y: 0
        });
      }
    }
    return {
      slug: id,
      area: true,
      key: `${center}/${width}${row.data("ftype")[1]}`,
      values: vals,
      peak: +center,
      type: row.data("ftype"),
      color: waveToColor(center),
      minwave: minX,
      maxwave: maxX
    };
  }

  function updateCustomFilter(row) {
    const F = customFilter(row);
    removeItem(row.attr("id"));
    data.push(F);
    refreshChart();
  }

  function updateCustomLaser(row) {
    const F = customLaser(row);
    removeItem(row.attr("id"));
    data.push(F);
    refreshChart();
  }

  // / Form Events

  $(".addFormItem").click(function onFormItemClick() {
    addFormItem(this.value, $(this).data("stype"), true);
  });

  let focusedItem;
  $("body").on("focus", ".data-selector, .select2", function onSelectFocus() {
    // Store the current value on focus and on change
    focusedItem = $(this)
      .closest(".row")
      .find(".data-selector")
      .val();
  });

  // main function when data-selector has been changed
  $("body").on("change", ".data-selector", function dataSelectorChange() {
    const selector = this;
    const slug = $(this).val();
    const row = $(this).closest(".row");
    // if the item is already selected, cancel operation
    if (
      dataHasSlug(slug) &&
      slug != focusedItem &&
      typeof localData[slug] !== "undefined"
    ) {
      alert(
        `${localData[slug][0].key
          .replace(" em", "")
          .replace(" 2p", "")
          .replace(" ex", "")} is already selected.`
      );
      if (focusedItem) {
        row
          .find(".data-selector")
          .val(focusedItem)
          .change();
      } else {
        row
          .find(".data-selector")
          .val(0)
          .change();
      }

      return;
    }

    // special behavior for custom bandpass
    if (this.value == "custom_bp") {
      row.children(":first").removeClass("col");
      row.children(":first").addClass("col-4");
      row.find(".custom_bp_form").show();
    } else {
      row.children(":first").removeClass("col-4");
      row.children(":first").addClass("col");
      row.find(".custom_bp_form").hide();
    }
    // special behavior for custom bandpass
    if (this.value == "custom_laser") {
      row.find(".custom_laser_form").show();
    } else {
      row.find(".custom_laser_form").hide();
    }

    // Remove the previous item fom the list
    removeItem(focusedItem);
    $(selector)
      .siblings(".item-link")
      .remove();
    // different process if it was a custom filter
    if (focusedItem == "custom_bp") {
      removeItem(row.attr("id"));
    }

    if (slug === "custom_bp") {
      updateCustomFilter(row);
    } else if (slug === "custom_laser") {
      updateCustomLaser(row);
    }
    // Add the new item to the data (unless it's blank)
    // then update the chart
    else if (slug !== null && slug !== 0) {
      addItem(slug).then(() => {
        if (localData[slug][0].url) {
          $(selector)
            .parent()
            .append(
              $("<div>", {
                class: "input-group-append item-link",
                title: "visit item page"
              }).append(
                $("<a>", {
                  class: "item-link",
                  href: localData[slug][0].url,
                  target: "_blank"
                }).append(
                  $("<button>", {
                    class: "btn btn-sm btn-info",
                    type: "button"
                  }).html($("#linkSVG").html())
                )
              )
            );
        }

        // if the slug has a two-2 item...
        if (localData[slug].map(e => e.type).indexOf(CONST.stype.twop) > -1) {
          // show the checkbox
          row.find("input.2pcheck").prop("disabled", false);
          // get the two photon item in the data
          const item2p = dataItemMatching({
            slug,
            type: CONST.stype.twop
          })[0]; // assume only 1
          // if the 2p button was already checked... keep 2p enabled
          if (row.find("input.2pcheck")[0].checked) {
            item2p.disabled = false;
          }
        } else {
          // if it doesn't have 2P data, hide the box...
          row.find("input.2pcheck").prop("disabled", true);
        }
        refreshChart();
      });
    } else {
      refreshChart();
    }

    focusedItem = slug;
  });

  $("body").on("click", ".remove-row", function onRemoveRow() {
    const row = $(this).closest(".row");
    let rowslug = row.find("select").val();
    if (rowslug === "custom_bp" || rowslug === "custom_laser") {
      rowslug = row.attr("id");
    }
    row.remove();
    removeItem(rowslug);

    refreshChart();
    if ($(".fluor-row").length <= 1) {
      $("#toggle_alls").hide();
    }
  });

  $("body").on("change", ".custom_bp_form input", e => {
    const $this = e.target;
    $this.value = Math.min($this.value, $this.max);
    $this.value = Math.max($this.value, $this.min);
    updateCustomFilter($($this).closest(".row"));
  });

  $("body").on("change", ".custom_laser_form input", e => {
    updateCustomLaser($(e.target).closest(".row"));
  });

  const chartsvg = $(selection)[0];
  const chartsvgNS = chartsvg.namespaceURI;
  const grad = document.createElementNS(chartsvgNS, "linearGradient");
  grad.setAttribute("id", "wavecolor_gradient");
  grad.setAttribute("class", "svggradient");
  const defs =
    chartsvg.querySelector("defs") ||
    chartsvg.insertBefore(
      document.createElementNS(chartsvgNS, "defs"),
      svg.firstChild
    );
  defs.appendChild(grad);

  function updateGlobalGradient() {
    $(grad).empty();
    const cmin = chart.xAxis.domain()[0];
    const cmax = chart.xAxis.domain()[1];
    const range = cmax - cmin;
    let stop;
    for (let w = cmin; w < cmax; w += 50) {
      stop = document.createElementNS(chartsvgNS, "stop");
      stop.setAttribute("offset", `${Math.round((100 * (w - cmin)) / range)}%`);
      stop.setAttribute("stop-color", waveToColor(w));
      grad.appendChild(stop);
    }
    stop = document.createElementNS(chartsvgNS, "stop");
    stop.setAttribute("offset", "100%");
    stop.setAttribute("stop-color", waveToColor(cmax));
    grad.appendChild(stop);
  }

  // // EFFICIENCY CALCULATIONS

  function spectraProduct(ar1, ar2) {
    // calculate product of two spectra.values
    const output = [];
    const left = Math.max(ar1[0].x, ar2[0].x);
    const right = Math.min(ar1[ar1.length - 1].x, ar2[ar2.length - 1].x);
    const offsetA1 = left - ar1[0].x; // these assume monotonic increase w/ step = 1
    const offsetA2 = left - ar2[0].x; // these assume monotonic increase w/ step = 1

    for (let i = 0; i < right - left; i++) {
      output.push({
        x: left + i,
        y: ar1[offsetA1 + i].y * ar2[offsetA2 + i].y
      });
    }
    return output;
  }

  function trapz(arr, min_, max_) {
    // approximate area under curve as series of trapezoids
    const min = min_ || 300;
    const max = max_ || 1000;
    let sum = 0;
    for (let i = 0; i < arr.length - 1; i++) {
      if (arr[i].x > min) {
        const d = (arr[i].y + arr[i + 1].y) / 2;
        sum += d;
      }
      if (arr[i].x > max) {
        break;
      }
    }
    return sum;
  }

  $("body").on("mousedown", ".effbutton", e => {
    const iSpec = $(e.target).data("ispec");
    const iFilt = $(e.target).data("ifilt");
    const overl = $(e.target).data("ioverlap");

    disabledData = [];
    for (let d = 0; d < data.length; d++) {
      disabledData.push(data[d].disabled);
      if (!(d === iFilt || d === iSpec)) {
        data[d].disabled = true;
      }
    }
    data.push({
      key: "Collection",
      values: overlaps[overl],
      color: "black",
      area: true
    });
    setTimeout(() => {
      chart.update();
    }, 60);
  });

  $("body").on("mouseup", ".effbutton", (e) => {
    const iSpec = $(e.target).data("ispec");
    const iFilt = $(e.target).data("ifilt");
    // const overl = $(e.target).data("ioverlap");

    data.splice(-1, 1);
    for (let d = 0; d < data.length; d++) {
      if (!((d === iFilt) || (d === iSpec))) {
        data[d].disabled = false;
      }
      data[d].disabled = disabledData[d];
    }

    //    chart.xDomain(getDomainLimits(data));
    refreshChart();
  });

  function eyebutton(emspect, emfilt, overlap) {
    return $("<button>", {
      class: "btn btn-info btn-sm mr-1 effbutton d-none d-sm-none d-md-inline",
      "data-ispec": emspect,
      "data-ifilt": emfilt,
      "data-ioverlap": overlap
    }).append($("#eyeSVG").html());
  }

  function calculateEfficiency() {
    const iSpectra = [];
    const iEmFilt = [];
    overlaps = [];
    $("#efficiency-table tbody").empty();
    $("#efficiency-table thead").empty();
    for (let i = 0; i < data.length; i++) {
      // look through current data for iSpectra and iEmFilt
      if (data[i].type == CONST.stype.em) {
        iSpectra.push(i);
      } else if (
        $.inArray(data[i].type, ["em_filter", "bp", "bm", "bs", "sp", "lp"]) >=
        0
      ) {
        iEmFilt.push(i);
        $("<th>")
          .text(data[i].key)
          .appendTo($("#efficiency-table").find("thead tr"));
      }
    }

    if (iSpectra.length && iEmFilt.length) {
      $("#efftab_blurb").hide();
      $("#efficiency-table thead").append($("<tr>").append($("<th>")));
      for (let x = 0; x < iEmFilt.length; x++) {
        $("<th>")
          .text(data[iEmFilt[x]].key)
          .appendTo($("#efficiency-table").find("thead tr"));
      }
      $(".efftab_help").show();
    } else {
      $("#efftab_blurb").show();
      $(".efftab_help").hide();
      return;
    }

    for (let s = 0; s < iSpectra.length; s++) {
      const emspectrum = data[iSpectra[s]];
      $(`<tr><td>${emspectrum.key}</td></tr>`).appendTo(
        $("#efficiency-table tbody")
      );
      for (let n = 0; n < iEmFilt.length; n++) {
        const EMtrans = spectraProduct(
          data[iEmFilt[n]].values,
          emspectrum.values
        );
        overlaps.push(EMtrans);
        const absEM = trapz(EMtrans);
        const EMpower = absEM / trapz(emspectrum.values);
        const formatted = Math.round(EMpower * 10000) / 100;
        let effclass = "efficiency-vbad";
        if (formatted > 75) {
          effclass = "efficiency-vgood";
        } else if (formatted > 50) {
          effclass = "efficiency-good";
        } else if (formatted > 25) {
          effclass = "efficiency-bad";
        }

        $("#efficiency-table tbody")
          .find("tr:last")
          .append(
            $("<td>", {
              class: effclass
            })
            .append(eyebutton(iSpectra[s], iEmFilt[n], overlaps.length - 1))
            .append(`  ${Math.round(absEM * 100) / 100} / (${formatted}%)`)
          );
      }
    }
  }

  // / IMORT MODALS

  $("#chromaImportForm, #semrockImportForm").submit(e => {
    e.preventDefault(); // avoid to execute the actual submit of the form.
    $("#footerSpinner").show();
    $("#footerFail").hide();
    $("#footerSuccess").hide();
    const form = $(e.target).closest("form");
    const brand = form.data("brand");
    $.ajax({
      type: "POST",
      url: form.attr("data-action-url"),
      data: form.serialize(),
      dataType: "json",
      success(data_) {
        if (data_.status) {
          const newdata = JSON.parse(data_.spectra_options);
          $('.data-selector[data-category="f"]').append(
            $("<option>", {
              value: newdata.slug
            }).text(newdata.name)
          );
          $(`#${brand}Input`).removeClass("is-invalid");
          $(`#${brand}Help`)
            .removeClass("invalid-feedback")
            .addClass("text-muted")
            .text("Success!");
          spectra_options.push(newdata); // override global variable with new options
          spectra_options.sort((a, b) =>
            a.name > b.name ? 1 : b.name > a.name ? -1 : 0
          );
          $("#footerSpinner").hide();
          $("#footerFail").hide();
          $("#footerSuccess").show();
        } else {
          $(`#${brand}Input`).addClass("is-invalid");
          $(`#${brand}Help`)
            .removeClass("text-muted")
            .addClass("invalid-feedback")
            .text(`ERROR: ${data_.message}`)
            .show();
          $("#footerSpinner").hide();
          $("#footerFail").show();
          $("#footerSuccess").hide();
        }
      }
    }).then(() => {
      $("#footerSpinner").hide();
      // $('#importModal').modal('hide');
    });
  });

  $(".importerClose").click(() => {
    $("#footerSpinner").hide();
    $("#footerFail").hide();
    $("#footerSuccess").hide();
  });

  // KEYBINDINGS

  function activateTab(tab) {
    $(`.nav-tabs a[href="#${tab}"]`).tab("show");
  }

  function findEmptyOrAdd(table, formtype, subtype) {
    focusedItem = 0;
    const el = $(`#${table}-table`).find('select option[value="0"]:selected');
    if (el.length) {
      el.last()
        .closest("select")
        .select2("open");
    } else {
      addFormItem(formtype, subtype, true);
    }
  }

  $("body").keypress(event => {
    // no double-events
    if ($(":focus").hasClass("select2-search__field")) {
      return;
    }

    if (document.activeElement.nodeName !== "INPUT") {
      // ALL KEYS REQUIRE SHIFT
      if (event.which === 80) {
        // p key
        activateTab("proteintab");
        findEmptyOrAdd("protein", "p");
      } else if (event.which === 68) {
        // d key
        activateTab("proteintab");
        findEmptyOrAdd("dye", "d");
      } else if (event.which === 77) {
        // m key
        activateTab("emtab");
        findEmptyOrAdd("emfilter", "f", "bm");
      } else if (event.which === 88) {
        // x key
        activateTab("extab");
        findEmptyOrAdd("exfilter", "f", "bx");
      } else if (event.which === 67) {
        // c key
        activateTab("emtab");
        findEmptyOrAdd("camqe", "c", "qe");
      } else if (event.which === 76) {
        // l key
        activateTab("extab");
        findEmptyOrAdd("light", "l");
      } else if (event.which === 79) {
        // o key
        activateTab("optionstab");
      } else if (event.which === 70) {
        // f key
        activateTab("efftab");
      }
    }
  });

  $("body").keydown(event => {
    if (event.which === 8) {
      // delete key
      const nextobject = $(":focus")
        .closest(".row")
        .prev(".row")
        .find(".select2-selection--single");
      $(":focus")
        .closest(".input-group")
        .find("button.remove-row")
        .click();
      nextobject.focus();
    }
    if (event.which === 38) {
      // up key
      event.preventDefault();
      $(":focus")
        .closest(".row")
        .prev(".row")
        .find(".select2-selection--single")
        .focus();
    }
    if (event.which === 40) {
      // up key
      event.preventDefault();
      $(":focus")
        .closest(".row")
        .next(".row")
        .find(".select2-selection--single")
        .focus();
    }

    if (event.which === 39) {
      // right key
      event.preventDefault();
      $(".nav-tabs .nav-link.active")
        .closest("li.nav-item")
        .next(".nav-item")
        .find(".nav-link")
        .click();
    }
    if (event.which === 37) {
      // left key
      event.preventDefault();
      $(".nav-tabs .nav-link.active")
        .closest("li.nav-item")
        .prev(".nav-item")
        .find(".nav-link")
        .click();
    }
  });

  return {
    getData() {
      return data;
    }
  };
}

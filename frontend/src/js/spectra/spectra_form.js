import { options, userOptions, CONST } from "./constants.js";

function uniqueID() {
  const id = Math.random()
    .toString(36)
    .substr(2, 9);
  return `_${id}`;
}

export function formSelection(filter, spectra_options) {
  const selWidget = $("<select>", {
    class: "custom-select custom-select-sm data-selector",
    "data-category": filter.category,
    "data-subtype": filter.subtype
  });
  selWidget.append(
    $("<option>", {
      value: 0
    }).text("-----")
  );
  if (filter.category == "f") {
    selWidget.append(
      $("<option>", {
        value: "custom_bp"
      }).text("Custom Bandpass")
    );
  }
  if (filter.category == "l") {
    selWidget.append(
      $("<option>", {
        value: "custom_laser"
      }).text("Laser Line")
    );
  }

  const matches = spectra_options.filter(function (item) {
    for (var key in filter) {
      if (typeof item[key] === 'undefined' || item[key] != filter[key]) {
        // FIXME
        // terrible hack to add bp filters to both excitation and emission dropdowns
        if ((filter[key] == 'bx' || filter[key] == 'bm') && item['subtype'] == 'bp') {
          return true;
        }
        if (filter[key] == 'bx' && item['subtype'] == 'sp') {
          return true;
        } if (filter[key] == 'bm' && item['subtype'] == 'lp') {
          return true;
        }
        return false;
      }
    }
    return true;
  });

  const usedSlugs = [];
  for (let i = 0; i < matches.length; i++) {
    if (usedSlugs.indexOf(matches[i].slug) == -1) {
      selWidget.append(
        $("<option>", {
          value: matches[i].slug
        }).text(matches[i].name)
      );
      usedSlugs.push(matches[i].slug);
    }
  }

  return selWidget;
}

export function fluorRow(widget) {
  const rowID = `f${  uniqueID()}`;
  return $("<div>", {
    class: "form-inline fluor-row row"
  })
    .attr("id", rowID)
    .append(
      $("<div>", {
        class: "col input-group input-group-sm"
      })
        .append(
          $("<div>", {
            class: "input-group-prepend select2-bootstrap-prepend"
          }).append(
            $("<button>", {
              class: "btn btn-danger btn-sm remove-row",
              type: "button"
            }).html("<strong>&times;</strong>")
          )
        )
        .append(widget.attr("id", `${rowID  }_selector`))
    )
    .append(
      $("<div>", {
        class: "col-1 custom-control custom-checkbox"
      })
        .append(
          $("<input>", {
            class: "excheck singlecheck custom-control-input",
            value: CONST.stype.ex,
            type: "checkbox",
            checked: "checked"
          }).attr("id", `${rowID  }_ex`)
        )
        .append(
          $("<label>", {
            class: "custom-control-label",
            for: `${rowID  }_ex`
          }).html('<span class="d-none d-sm-inline">e</span>x')
        )
    )
    .append(
      $("<div>", {
        class: "col-1 custom-control custom-checkbox"
      })
        .append(
          $("<input>", {
            class: "emcheck singlecheck custom-control-input",
            value: CONST.stype.em,
            type: "checkbox",
            checked: "checked"
          }).attr("id", `${rowID  }_em`)
        )
        .append(
          $("<label>", {
            class: "custom-control-label",
            for: `${rowID  }_em`
          }).html('<span class="d-none d-sm-inline">e</span>m')
        )
    )
    .append(
      $("<div>", {
        class: "col-1 custom-control custom-checkbox"
      })
        .append(
          $("<input>", {
            class: "2pcheck singlecheck custom-control-input",
            value: CONST.stype.twop,
            type: "checkbox",
            checked: "checked"
          })
            .attr("id", `${rowID  }_2p`)
            .prop("checked", !options.hide2p)
        )
        .append(
          $("<label>", {
            class: "custom-control-label",
            for: `${rowID  }_2p`
          }).html('2<span class="d-none d-sm-inline">p</span>')
        )
    );
}

export function excRow(widget, cls, refreshChart) {
  const rowID = `l${  uniqueID()}`;
  return $("<div>", {
    class: `form-inline row ${  cls}`
  })
    .attr("id", rowID)
    .append(
      $("<div>", {
        class: "col input-group input-group-sm"
      })
        .append(
          $("<div>", {
            class: "input-group-prepend select2-bootstrap-prepend"
          }).append(
            $("<button>", {
              class: "btn btn-danger btn-sm remove-row",
              type: "button"
            })
              .html("<strong>&times;</strong>")
              .click(function() {
                options.exNormWave = undefined;
              })
          )
        )
        .append(
          widget.change(function() {
            let normchecked = $(this)
              .closest(`.${  cls}`)
              .find(".exnormcheck")
              .prop("checked");
            let islaser = $(this).val() == "custom_laser";
            if (islaser) {
              if (normchecked) {
                options.exNormWave = +$(this)
                  .parent()
                  .next()
                  .find(".custom_laser_wave")
                  .val();
                refreshChart();
              }
            } else {
              options.exNormWave = undefined;
              refreshChart();
            }
          })
        )
    )
    .append(
      $("<div>", {
        class: "col-3 input-group input-group-sm custom_laser_form"
      })
        .hide()
        .append(
          $("<div>", {
            class: "input-group-prepend"
          }).append(
            $("<span>", {
              class: "input-group-text",
              id: `${rowID  }_laserwave`
            }).html("&lambda;")
          )
        )
        .append(
          $("<input>", {
            class: "form-control form-control-sm custom_laser_wave",
            type: "number",
            name: "custom_laser_wave",
            min: "300",
            max: "1500",
            value: "488",
            "aria-label": "Laser Wavelength",
            "aria-describedby": `${rowID  }_laserwave`
          }).change(function() {
            if (
              $(this)
                .closest(".row")
                .find(".exnormcheck")
                .prop("checked")
            ) {
              options.exNormWave = +this.value;
              refreshChart();
            }
          })
        )
    )
    .append(
      $("<div>", {
        class: "col-2 custom-control custom-radio ml-2 custom_laser_form"
      })
        .hide()
        .append(
          $("<input>", {
            class: "exnormcheck custom-control-input",
            "data-checktype": "exnorm",
            "data-lastval": " ",
            type: "radio",
            name: "exnormRadio",
            value: rowID,
            id: `${rowID  }_norm`
          }).click(function(e) {
            let lastval = $(this).data("lastval");
            $(".exnormcheck").data("lastval", "");
            if (this.value == lastval) {
              // already clicked, unlick
              $("#exnormRadioOFF").prop("checked", true);
              options.exNormWave = undefined;
            } else {
              options.exNormWave = undefined;
            const v = $(this).closest(`.${  cls}`).find('.data-selector').val();
              if (v == "custom_laser") {
                options.exNormWave = +$(this)
                  .parent()
                  .prev()
                  .find(".custom_laser_wave")
                  .val();
              }
              $(this).data("lastval", this.value);
            }
            refreshChart();
          })
        )
        .append(
          $("<label>", {
            class: "custom-control-label",
            for: `${rowID  }_norm`
          }).text("norm em")
        )
    );
}

export function filterRow(widget, stype) {
  const rowID = stype + uniqueID();
  return $("<div>", {
    class: "form-inline filter-row row",
    "data-ftype": stype
  })
    .attr("id", rowID)
    .append(
      $("<div>", {
        class: "col-6 col-md input-group input-group-sm"
      })
        .append(
          $("<div>", {
            class: "input-group-prepend select2-bootstrap-prepend"
          }).append(
            $("<button>", {
              class: "btn btn-danger btn-sm remove-row",
              type: "button"
            }).html("<strong>&times;</strong>")
          )
        )
        .append(widget.attr("id", `${rowID  }_selector`))
    )
    .append(
      $("<div>", {
        class: "col-6 col-md-3 input-group input-group-sm custom_bp_form"
      })
        .hide()
        .append(
          $("<div>", {
            class: "input-group-prepend"
          }).append(
            $("<span>", {
              class: "input-group-text",
              id: `${rowID  }_centerwave`
            }).text("center")
          )
        )
        .append(
          $("<input>", {
            class: "form-control form-control-sm custom_em_value",
            type: "number",
            name: "custom_em_center",
            min: "300",
            max: "1500",
            value: "525",
            "aria-label": "Filter Center",
            "aria-describedby": `${rowID  }_centerwave`
          })
        )
    )
    .append(
      $("<div>", {
        class: "col-6 col-md-3 input-group input-group-sm custom_bp_form"
      })
        .hide()
        .append(
          $("<div>", {
            class: "input-group-prepend"
          }).append(
            $("<span>", {
              class: "input-group-text",
              id: `${rowID  }_bandwidth`
            }).text("width")
          )
        )
        .append(
          $("<input>", {
            class: "form-control form-control-sm custom_em_value",
            type: "number",
            name: "custom_em_bandwidth",
            min: "1",
            max: "800",
            value: "50",
            "aria-label": "Filter Width",
            "aria-describedby": `${rowID  }_bandwidth`
          })
        )
    )
    .append(
      $("<div>", {
        class: "col-6 col-md-2 input-group input-group-sm custom_bp_form"
      })
        .hide()
        .append(
          $("<div>", {
            class: "input-group-prepend"
          }).append(
            $("<span>", {
              class: "input-group-text",
              id: `${rowID  }_trans`
            }).text("%T")
          )
        )
        .append(
          $("<input>", {
            class: "form-control form-control-sm custom_em_value",
            type: "number",
            name: "custom_em_trans",
            min: "0.001",
            max: "1",
            value: "0.95",
            step: "0.01",
            "aria-label": "Filter Width",
            "aria-describedby": `${rowID  }_trans`
          })
        )
    );
}

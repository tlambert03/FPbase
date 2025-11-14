const $ = window.jQuery // jQuery loaded from CDN

export default function initSearch(filterfields, operatorLookup, labelLookup) {
  var fields = {}

  for (var key in filterfields) {
    if (Object.hasOwn(filterfields, key)) {
      fields[key] = new Set(filterfields[key])
    }
  }

  var available_fields = new Set(Object.keys(fields))
  var i = 0
  var prevName
  var prevOperator

  function toTitleCase(str) {
    if (!str) {
      return ""
    }
    str = str.split("__").join("_")
    str = str.split("_").join(" ")
    return str.replace(/\w\S*/g, (txt) => txt.charAt(0).toUpperCase() + txt.substr(1).toLowerCase())
  }

  function name_options() {
    var out = ""
    for (const n of available_fields) {
      if (n !== "id" && n !== "slug") {
        let name
        if (n in labelLookup) {
          name = labelLookup[n]
        } else {
          name = toTitleCase(n)
        }
        out += `<option value=${n}>${name}</option>`
      }
    }
    return out
  }

  function operator_options(ops) {
    var out = ""
    for (const o of ops) {
      let op
      if (o in operatorLookup) {
        op = operatorLookup[o]
      } else {
        op = toTitleCase(ops[o])
      }
      out += `<option value=${o}>${op}</option>`
    }
    return out
  }

  $.fn.exists = function () {
    return this.length !== 0
  }

  function updateInputField(row, keepValue) {
    // figure out filter name and operator
    var filter_name = row.find(".filter-select").val()
    var operation = row.find(".operator-select").val()
    let name
    if (operation === "exact") {
      name = filter_name
    } else {
      name = `${filter_name}__${operation}`
    }

    // hide the current input field by putting it in the hidden #crispy-form
    // FIXME: change this to a variable selector
    var inputcol = row.find(".input-col")
    let value
    if (keepValue) {
      value = inputcol.find("input").val()
    }
    inputcol.find(".form-group").appendTo("#crispy-form")

    // grab the form div that we want to put here
    var formdiv = $(`#div_id_${name}`)
    formdiv.appendTo(inputcol)

    // for some reason range input is losing this class
    formdiv.find("input").addClass("form-control")
    //formdiv.find('input').prop('required',true);
    if (keepValue) {
      formdiv.find("input").val(value)
    }
  }

  function queryRow(id) {
    var out =
      '\
      <div class="row g-2 query-row" id="query-row-' +
      id +
      '">\
          <div class="col-sm-4 names-col d-flex justify-content-between align-items-start">\
          <button type="button" class="btn btn-danger remove-row-btn form-group me-2" ><strong>&times;</strong></button>\
                  <div class="form-group" style="width:100%;"><select class="form-control filter-select" id="filter-select-' +
      id +
      '">' +
      name_options() +
      '</select>\
              </div>\
          </div>\
          <div class="col-sm-4 operator-col">\
              <div class="form-group"><select class="form-control operator-select">\
              </select></div>\
          </div>\
          <div class="col-sm-4">\
                  <div class="form-group input-col"></div>\
          </div>\
      </div>'
    return out
  }

  function addRow(target, filter, operator) {
    const newrow = $(queryRow(i)).appendTo(target)

    if (filter) {
      newrow.find(".filter-select").val(filter)
    }

    const filterSelector = newrow.find(".filter-select")
    const filterName = filterSelector.val()

    const operatorSelect = newrow.find(".operator-select")
    operatorSelect.html(operator_options(fields[filterName]))
    operatorSelect.removeClass()
    operatorSelect.addClass(`form-control operator-select ${filterName}_operator`)

    if (operator) {
      operatorSelect.val(operator)
    } else {
      operator = operatorSelect.val()
    }

    disableOperator(filterName, operator)
    updateOperators(filterName, operatorSelect)
    updateInputField(newrow)

    //updateOperators(row.find('.filter-select').val());
    i += 1
  }

  function disableOperator(name, op) {
    fields[name].delete(op)
    if (fields[name].size === 0) {
      available_fields.delete(name)
    }
  }

  function enableOperator(name, op) {
    fields[name].add(op)
    available_fields.add(name)
  }

  function updateOperators(filterName, sender) {
    // where sender is the operatorSelect that sent the update command
    const selector = $(`.${filterName}_operator`).not(sender)
    if (selector.length > 0) {
      selector.find("option").not(":selected").remove()
      selector.append(operator_options(fields[filterName]))
    }
  }

  $("body").on("focus", ".filter-select", function (_event) {
    // Store the current value on focus and on change
    prevName = this.value
    prevOperator = $(this).closest(".query-row").find(".operator-select").val()
  })
  $("body").on("focus", ".operator-select", function (_event) {
    // Store the current value on focus and on change
    prevOperator = this.value
  })

  $("#add-row-btn").on("click", () => {
    addRow("#query_builder")
  })

  $("body").on("click", ".remove-row-btn", function () {
    const thisrow = $(this).closest(".query-row")
    const filterName = thisrow.find(".filter-select").val()
    const operatorSelect = thisrow.find(".operator-select")
    enableOperator(filterName, operatorSelect.val())
    updateOperators(filterName, operatorSelect)

    const inputfield = thisrow.find(".input-col").find(".form-group")
    inputfield.appendTo("#crispy-form")

    thisrow.remove()
    i -= 1
  })

  $("body").on("change", ".filter-select", function () {
    //when the filter name selector gets changed (or added)
    const dropdownSelect = $(this)
    const filterName = dropdownSelect.val()

    const thisrow = $(this).closest(".query-row")

    // remove the old operator dropdown and add the new one
    const operatorSelect = thisrow.find(".operator-select")
    operatorSelect.empty()
    operatorSelect.html(operator_options(fields[filterName]))
    operatorSelect.removeClass()
    operatorSelect.addClass(`form-control operator-select ${filterName}_operator`)

    if (prevName && prevOperator) {
      enableOperator(prevName, prevOperator)
    }
    disableOperator(filterName, operatorSelect.val())
    updateOperators(prevName)
    updateOperators(filterName, operatorSelect)
    updateInputField(thisrow)

    prevName = this.value
    prevOperator = $(this).closest(".query-row").find(".operator-select").val()
  })

  $("body").on("change", ".operator-select", function () {
    const thisrow = $(this).closest(".query-row")
    const filterName = thisrow.find(".filter-select").val()
    enableOperator(filterName, prevOperator)
    disableOperator(filterName, this.value)
    updateOperators(filterName)
    updateInputField(thisrow, true)

    prevOperator = this.value
  })

  function loadState(state) {
    for (const key in state) {
      const value = state[key]
      if (key === "display") {
        $(`#${value}button`).click()
        continue
      }
      let filter, operator
      if (key in fields) {
        filter = key
        operator = "exact"
      } else {
        const splits = key.split("__")
        filter = splits.slice(0, splits.length - 1).join("__")
        operator = splits[splits.length - 1]
      }
      if (filter && operator) {
        addRow("#query_builder", filter, operator)
      }
    }
  }

  $(() => {
    $("#crispy-form label").remove()

    var urlParams
    window.onpopstate = () => {
      var match,
        pl = /\+/g, // Regex for replacing addition symbol with a space
        search = /([^&=]+)=?([^&]*)/g,
        decode = (s) => decodeURIComponent(s.replace(pl, " ")),
        query = window.location.search.substring(1)

      urlParams = {}
      match = search.exec(query)
      while (match) {
        urlParams[decode(match[1])] = decode(match[2])
        match = search.exec(query)
      }
    }
    window.onpopstate()

    $("#query_builder").empty()
    loadState(urlParams)
    if (!$("#query_builder").find(".query-row").exists()) {
      $("#query_builder").empty()
      addRow("#query_builder")
    }

    $(".displaybuttons input").change(function () {
      const display_type = $(this).val()
      $(`#${display_type}display`).removeClass("hidden")
      $(`#${display_type}display`).siblings("div").addClass("hidden")
    })

    // Mark the form as ready for E2E tests
    $("#query_builder").attr("data-search-ready", "true")
  })
}

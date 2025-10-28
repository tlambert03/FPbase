window.initSearch = function(filterfields, operatorLookup, labelLookup) {
  var fields = {};

  for (var key in filterfields) {
      if (filterfields.hasOwnProperty(key)) {
        fields[key] = new Set(filterfields[key]);
      }
  }

  var available_fields = new Set(Object.keys(fields));
  var i = 0;
  var prevName;
  var prevOperator;

  function toTitleCase(str) {
    if (!str) {
      return '';
    }
    str = str.split('__').join('_');
    str = str.split('_').join(' ');
    return str.replace(/\w\S*/g, function(txt) {
      return txt.charAt(0).toUpperCase() + txt.substr(1).toLowerCase();
    });
  }

  function name_options() {
    var out = '';
    for (let n of available_fields) {
      if (n !== 'id' && n !== 'slug') {
        if (n in labelLookup) {
          name = labelLookup[n];
        } else {
          name = toTitleCase(n);
        }
        out += '<option value=' + n + '>' + name + '</option>';
      }
    }
    return out;
  }

  function operator_options(ops) {
    var out = '';
    for (let o of ops) {
      if (o in operatorLookup) {
        op = operatorLookup[o];
      } else {
        op = toTitleCase(ops[o]);
      }
      out += '<option value=' + o + '>' + op + '</option>';
    }
    return out;
  }

  $.fn.exists = function() {
    return this.length !== 0;
  };

  function updateInputField(row, keepValue) {
    // figure out filter name and operator
    var filter_name = row.find('.filter-select').val();
    var operation = row.find('.operator-select').val();
    if (operation == 'exact') {
      name = filter_name;
    } else {
      name = filter_name + '__' + operation;
    }

    // hide the current input field by putting it in the hidden #crispy-form
    // FIXME: change this to a variable selector
    var inputcol = row.find('.input-col');
    if (keepValue) {
      value = inputcol.find('input').val();
    }
    inputcol.find('.form-group').appendTo('#crispy-form');

    // grab the form div that we want to put here
    var formdiv = $('#div_id_' + name);
    formdiv.appendTo(inputcol);

    // for some reason range input is losing this class
    formdiv.find('input').addClass('form-control');
    //formdiv.find('input').prop('required',true);
    if (keepValue) {
      formdiv.find('input').val(value);
    }
  }

  function queryRow(id) {
    var out =
      '\
      <div class="form-row query-row" id="query-row-' +
      id +
      '">\
          <div class="col-sm-4 names-col d-flex justify-content-between align-items-start">\
          <button type="button" class="btn btn-danger remove-row-btn form-group mr-2" ><strong>&times;</strong></button>\
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
      </div>';
    return out;
  }

  function addRow(target, filter, operator) {
    newrow = $(queryRow(i)).appendTo(target);

    if (filter) {
      newrow.find('.filter-select').val(filter);
    }

    filterSelector = newrow.find('.filter-select');
    filterName = filterSelector.val();

    operatorSelect = newrow.find('.operator-select');
    operatorSelect.html(operator_options(fields[filterName]));
    operatorSelect.removeClass();
    operatorSelect.addClass(
      'form-control operator-select ' + filterName + '_operator'
    );

    if (operator) {
      operatorSelect.val(operator);
    } else {
      operator = operatorSelect.val();
    }

    disableOperator(filterName, operator);
    updateOperators(filterName, operatorSelect);
    updateInputField(newrow);

    //updateOperators(row.find('.filter-select').val());
    i += 1;
  }

  function disableOperator(name, op) {
    fields[name].delete(op);
    if (fields[name].size == 0) {
      available_fields.delete(name);
    }
  }

  function enableOperator(name, op) {
    fields[name].add(op);
    available_fields.add(name);
  }

  function updateOperators(filterName, sender) {
    // where sender is the operatorSelect that sent the update command
    selector = $('.' + filterName + '_operator').not(sender);
    if (selector.length > 0) {
      selector
        .find('option')
        .not(':selected')
        .remove();
      selector.append(operator_options(fields[filterName]));
    }
  }

  $('body').on('focus', '.filter-select', function(event) {
    // Store the current value on focus and on change
    prevName = this.value;
    prevOperator = $(this)
      .closest('.query-row')
      .find('.operator-select')
      .val();
  });
  $('body').on('focus', '.operator-select', function(event) {
    // Store the current value on focus and on change
    prevOperator = this.value;
  });

  $('#add-row-btn').on('click', function() {
    addRow('#query_builder');
  });

  $('body').on('click', '.remove-row-btn', function() {
    thisrow = $(this).closest('.query-row');
    filterName = thisrow.find('.filter-select').val();
    operatorSelect = thisrow.find('.operator-select');
    enableOperator(filterName, operatorSelect.val());
    updateOperators(filterName, operatorSelect);

    inputfield = thisrow.find('.input-col').find('.form-group');
    inputfield.appendTo('#crispy-form');

    thisrow.remove();
    i -= 1;
  });

  $('body').on('change', '.filter-select', function() {
    //when the filter name selector gets changed (or added)
    dropdownSelect = $(this);
    filterName = dropdownSelect.val();

    thisrow = $(this).closest('.query-row');

    // remove the old operator dropdown and add the new one
    operatorSelect = thisrow.find('.operator-select');
    operatorSelect.empty();
    operatorSelect.html(operator_options(fields[filterName]));
    operatorSelect.removeClass();
    operatorSelect.addClass(
      'form-control operator-select ' + filterName + '_operator'
    );

    if (prevName && prevOperator) {
      enableOperator(prevName, prevOperator);
    }
    disableOperator(filterName, operatorSelect.val());
    updateOperators(prevName);
    updateOperators(filterName, operatorSelect);
    updateInputField(thisrow);

    prevName = this.value;
    prevOperator = $(this)
      .closest('.query-row')
      .find('.operator-select')
      .val();
  });

  $('body').on('change', '.operator-select', function() {
    thisrow = $(this).closest('.query-row');
    filterName = thisrow.find('.filter-select').val();
    enableOperator(filterName, prevOperator);
    disableOperator(filterName, this.value);
    updateOperators(filterName);
    updateInputField(thisrow, true);

    prevOperator = this.value;
  });

  function loadState(state) {
    for (key in state) {
      value = state[key];
      if (key === 'display') {
        $('#' + value + 'button').click();
        continue;
      }
      if (key in fields) {
        filter = key;
        operator = 'exact';
      } else {
        splits = key.split('__');
        filter = splits.slice(0, splits.length - 1).join('__');
        operator = splits[splits.length - 1];
      }
      if (filter && operator){
        addRow('#query_builder', filter, operator);
      }
    }
  }

  $(function() {
    $('#crispy-form label').remove();

    var urlParams;
    (window.onpopstate = function() {
      var match,
        pl = /\+/g, // Regex for replacing addition symbol with a space
        search = /([^&=]+)=?([^&]*)/g,
        decode = function(s) {
          return decodeURIComponent(s.replace(pl, ' '));
        },
        query = window.location.search.substring(1);

      urlParams = {};
      while ((match = search.exec(query)))
        urlParams[decode(match[1])] = decode(match[2]);
    })();

    $('#query_builder').empty();
    loadState(urlParams);
    if (
      !$('#query_builder')
        .find('.query-row')
        .exists()
    ) {
      $('#query_builder').empty();
      addRow('#query_builder');
    }

    $('.displaybuttons input').change(function() {
      var display_type = $(this).val();
      $('#' + display_type + 'display').show();
      $('#' + display_type + 'display')
        .siblings('div')
        .hide();
    });
  });
};

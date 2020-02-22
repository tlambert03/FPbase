/**
 * jQuery Formset 1.3-pre
 * @author Stanislaus Madueke (stan DOT madueke AT gmail DOT com)
 * @requires jQuery 1.2.6 or later
 *
 * Copyright (c) 2009, Stanislaus Madueke
 * All rights reserved.
 *
 * Licensed under the New BSD License
 * See: http://www.opensource.org/licenses/bsd-license.php
 */
;(function($) {
  $.fn.formset = function(opts) {
    const options = $.extend({}, $.fn.formset.defaults, opts)
    const flatExtraClasses = options.extraClasses.join(" ")
    const totalForms = $(`#id_${options.prefix}-TOTAL_FORMS`)
    const maxForms = $(`#id_${options.prefix}-MAX_NUM_FORMS`)
    const minForms = $(`#id_${options.prefix}-MIN_NUM_FORMS`)
    const childElementSelector = "input,select,textarea,label,div"
    const $$ = $(this)
    const applyExtraClasses = function(row, ndx) {
      if (options.extraClasses) {
        row.removeClass(flatExtraClasses)
        row.addClass(options.extraClasses[ndx % options.extraClasses.length])
      }
    }
    const updateElementIndex = function(elem, prefix, ndx) {
      const idRegex = new RegExp(`${prefix}-(\\d+|__prefix__)-`)
      const replacement = `${prefix}-${ndx}-`
      if (elem.attr("for"))
        elem.attr("for", elem.attr("for").replace(idRegex, replacement))
      if (elem.attr("id"))
        elem.attr("id", elem.attr("id").replace(idRegex, replacement))
      if (elem.attr("name"))
        elem.attr("name", elem.attr("name").replace(idRegex, replacement))
    }
    const hasChildElements = function(row) {
      return row.find(childElementSelector).length > 0
    }
    const showAddButton = function() {
      return (
        maxForms.length == 0 || // For Django versions pre 1.2
        maxForms.val() == "" ||
        maxForms.val() - totalForms.val() > 0
      )
    }
    /**
     * Indicates whether delete link(s) can be displayed - when total forms > min forms
     */
    const showDeleteLinks = function() {
      return (
        minForms.length == 0 || // For Django versions pre 1.7
        minForms.val() == "" ||
        totalForms.val() - minForms.val() > 0
      )
    }
    const insertDeleteLink = function(row) {
      const delCssSelector = $.trim(options.deleteCssClass).replace(/\s+/g, ".")
      const addCssSelector = $.trim(options.addCssClass).replace(/\s+/g, ".")
      if (row.is("TR")) {
        // If the forms are laid out in table rows, insert
        // the remove button into the last table cell:
        row
          .children(":last")
          .append(
            `<a class="${options.deleteCssClass}" href="javascript:void(0)">${options.deleteText}</a>`
          )
      } else if (row.is("UL") || row.is("OL")) {
        // If they're laid out as an ordered/unordered list,
        // insert an <li> after the last list item:
        row.append(
          `<li><a class="${options.deleteCssClass}" href="javascript:void(0)">${options.deleteText}</a></li>`
        )
      } else {
        // Otherwise, just insert the remove button as the
        // last child element of the form's container:
        row.append(
          `<a class="${options.deleteCssClass}" href="javascript:void(0)">${options.deleteText}</a>`
        )
      }
      // Check if we're under the minimum number of forms - not to display delete link at rendering
      if (!showDeleteLinks()) {
        row.find(`a.${delCssSelector}`).hide()
      }

      row.find(`a.${delCssSelector}`).click(function() {
        const row = $(this).parents(`.${options.formCssClass}`)
        const del = row.find('input:hidden[id $= "-DELETE"]')
        const buttonRow = row.siblings(
          `a.${addCssSelector}, .${options.formCssClass}-add`
        )
        let forms
        if (del.length) {
          // We're dealing with an inline formset.
          // Rather than remove this form from the DOM, we'll mark it as deleted
          // and hide it, then let Django handle the deleting:
          del.val("on")
          row.hide()
          forms = $(`.${options.formCssClass}`).not(":hidden")
        } else {
          row.remove()
          // Update the TOTAL_FORMS count:
          forms = $(`.${options.formCssClass}`).not(".formset-custom-template")
          totalForms.val(forms.length).trigger("change")
        }
        for (var i = 0, formCount = forms.length; i < formCount; i++) {
          // Apply `extraClasses` to form rows so they're nicely alternating:
          applyExtraClasses(forms.eq(i), i)
          if (!del.length) {
            // Also update names and IDs for all child controls (if this isn't
            // a delete-able inline formset) so they remain in sequence:
            forms
              .eq(i)
              .find(childElementSelector)
              .each(function() {
                updateElementIndex($(this), options.prefix, i)
              })
          }
        }
        // Check if we've reached the minimum number of forms - hide all delete link(s)
        if (!showDeleteLinks()) {
          $(`a.${delCssSelector}`).each(function() {
            $(this).hide()
          })
        }
        // Check if we need to show the add button:
        if (buttonRow.is(":hidden") && showAddButton()) buttonRow.show()
        // If a post-delete callback was provided, call it with the deleted form:
        if (options.removed) options.removed(row)
        return false
      })
    }

    $$.each(function(i) {
      const row = $(this)
      const del = row.find('input:checkbox[id $= "-DELETE"]')
      if (del.length) {
        // If you specify "can_delete = True" when creating an inline formset,
        // Django adds a checkbox to each form in the formset.
        // Replace the default checkbox with a hidden field:
        if (del.is(":checked")) {
          // If an inline formset containing deleted forms fails validation, make sure
          // we keep the forms hidden (thanks for the bug report and suggested fix Mike)
          del.before(
            `<input type="hidden" name="${del.attr("name")}" id="${del.attr(
              "id"
            )}" value="on" />`
          )
          row.hide()
        } else {
          del.before(
            `<input type="hidden" name="${del.attr("name")}" id="${del.attr(
              "id"
            )}" />`
          )
        }
        // Hide any labels associated with the DELETE checkbox:
        $(`label[for="${del.attr("id")}"]`).hide()
        del.remove()
      }
      if (hasChildElements(row)) {
        row.addClass(options.formCssClass)
        if (row.is(":visible") || options.processHidden) {
          insertDeleteLink(row)
          applyExtraClasses(row, i)
        }
      }
    })

    if ($$.length) {
      const hideAddButton = !showAddButton()
      let addButton
      let template
      if (options.formTemplate) {
        // If a form template was specified, we'll clone it to generate new form instances:
        template =
          options.formTemplate instanceof $
            ? options.formTemplate
            : $(options.formTemplate)
        template
          .removeAttr("id")
          .addClass(`${options.formCssClass} formset-custom-template`)
        template.find(childElementSelector).each(function() {
          updateElementIndex($(this), options.prefix, "__prefix__")
        })
        insertDeleteLink(template)
      } else {
        // Otherwise, use the last form in the formset; this works much better if you've got
        // extra (>= 1) forms (thnaks to justhamade for pointing this out):
        template = $(`.${options.formCssClass}:last`)
          .clone(true)
          .removeAttr("id")
        template.find('input:hidden[id $= "-DELETE"]').remove()
        // Clear all cloned fields, except those the user wants to keep (thanks to brunogola for the suggestion):
        template
          .find(childElementSelector)
          .not(options.keepFieldValues)
          .each(function() {
            const elem = $(this)
            // If this is a checkbox or radiobutton, uncheck it.
            // This fixes Issue 1, reported by Wilson.Andrew.J:
            if (elem.is("input:checkbox") || elem.is("input:radio")) {
              elem.attr("checked", false)
            } else {
              elem.val("")
            }
          })
      }
      // FIXME: Perhaps using $.data would be a better idea?
      options.formTemplate = template

      if ($$.is("TR")) {
        // If forms are laid out as table rows, insert the
        // "add" button in a new table row:
        const numCols = $$.eq(0).children().length // This is a bit of an assumption :|
        const buttonRow = $(
          `<tr><td colspan="${numCols}"><a class="${options.addCssClass}" href="javascript:void(0)">${options.addText}</a></tr>`
        ).addClass(`${options.formCssClass}-add`)
        $$.parent().append(buttonRow)
        if (hideAddButton) buttonRow.hide()
        addButton = buttonRow.find("a")
      } else {
        // Otherwise, insert it immediately after the last form:
        $$.filter(":last").after(
          `<a class="${options.addCssClass}" href="javascript:void(0)">${options.addText}</a>`
        )
        addButton = $$.filter(":last").next()
        if (hideAddButton) addButton.hide()
      }
      addButton.click(function() {
        const formCount = parseInt(totalForms.val())
        const row = options.formTemplate
          .clone(true)
          .removeClass("formset-custom-template")
        const buttonRow = $(
          $(this)
            .parents(`tr.${options.formCssClass}-add`)
            .get(0) || this
        )
        const delCssSelector = $.trim(options.deleteCssClass).replace(
          /\s+/g,
          "."
        )
        applyExtraClasses(row, formCount)
        row.insertBefore(buttonRow).show()
        row.find(childElementSelector).each(function() {
          updateElementIndex($(this), options.prefix, formCount)
        })
        totalForms.val(formCount + 1).trigger("change")
        // Check if we're above the minimum allowed number of forms -> show all delete link(s)
        if (showDeleteLinks()) {
          $(`a.${delCssSelector}`).each(function() {
            $(this).show()
          })
        }
        // Check if we've exceeded the maximum allowed number of forms:
        if (!showAddButton()) buttonRow.hide()
        // If a post-add callback was supplied, call it with the added form:
        if (options.added) options.added(row)
        return false
      })
    }

    return $$
  }

  /* Setup plugin defaults */
  $.fn.formset.defaults = {
    prefix: "form", // The form prefix for your django formset
    formTemplate: null, // The jQuery selection cloned to generate new form instances
    addText: "add another", // Text for the add link
    deleteText: "remove", // Text for the delete link
    addCssClass: "add-row", // CSS class applied to the add link
    deleteCssClass: "delete-row", // CSS class applied to the delete link
    formCssClass: "dynamic-form", // CSS class applied to each form in a formset
    extraClasses: [], // Additional CSS classes, which will be applied to each form in turn
    keepFieldValues: "", // jQuery selector for fields whose values should be kept when the form is cloned
    added: null, // Function called each time a new form is added
    removed: null, // Function called each time a form is deleted
  }
})(jQuery)

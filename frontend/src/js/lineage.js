import d3 from "d3"
import $ from "jquery"

function n_sibs(node) {
  return node.parent ? node.parent.children.length - 1 : 0
}

function is_first_sib(node) {
  if (!node.parent || node.parent.children.length === 1) {
    return false
  }
  var i = node.parent.children
    .map(function(e) {
      return e.name
    })
    .indexOf(node.name)
  if (i === 0) {
    return true
  }
  return false
}

function is_last_sib(node) {
  if (!node.parent || node.parent.children.length === 1) {
    return false
  }
  var i = node.parent.children
    .map(function(e) {
      return e.name
    })
    .indexOf(node.name)
  if (i === node.parent.children.length - 1) {
    return true
  }
  return false
}

function text_position(node, slug, vertical) {
  var x
  var y
  var anchor = "middle"
  if (node.depth < 2) {
    x = -13
    y = 0
    anchor = "end"
  } else if (n_sibs(node) < 2) {
    if (node.children || node._children) {
      x = 0
      y = is_last_sib(node) ? 16 : -16
    } else {
      x = 13
      y = 0
      anchor = "start"
    }
  } else {
    y = 0
    if (node.children || node._children) {
      if (is_first_sib(node)) {
        y = -16
        x = 0
        anchor = "middle"
      } else if (is_last_sib(node)) {
        y = 16
        x = 0
        anchor = "middle"
      } else {
        x = -13
        anchor = "end"
      }
    } else {
      x = 13
      anchor = "start"
    }
  }
  if (node.slug === slug) {
    if (node.children || node._children) {
      y -= 8
    } else {
      x += 6
    }
  }
  if (vertical) {
    return [y, x, anchor]
  }
  return [x, y, anchor]
}

export default function LineageChart(conf) {
  let config = conf || {}
  let margin = config.margin || { top: 20, right: 110, bottom: 15, left: 65 },
    width,
    minNodeWidth = config.minNodeWidth || 70,
    heightScalar = 40,
    widthScalar = 1,
    nodeWidth,
    height,
    i = 0,
    duration = config.duration || 0,
    defaultRadius = 8,
    slugRadius = 13,
    data,
    sel,
    slug = config.slug || null,
    show_inserts = true,
    show_deletions = true,
    tree = d3.layout.tree(),
    withTopScroll = config.withTopScroll || false,
    withSearch = config.withSearch || false,
    withToolbar = config.withToolbar || false,
    vertical = config.vertical || false,
    dropShadow = {
      stdDeviation: 4,
      dx: 0,
      dy: 0,
      floodColor: "#38e"
    }

  var diagonal = d3.svg.diagonal().projection(function(d) {
    return [d.y, d.x]
  })

  // Define the div for the tooltip
  var tooltip = d3
    .select("body")
    .append("div")
    .attr("class", "tooltip lineage-tooltip")
    .style("opacity", 0)

  // function stopEvent(event) {
  //   if (event.preventDefault !== undefined) event.preventDefault();
  //   if (event.stopPropagation !== undefined) event.stopPropagation();
  // }

  function chart(selection) {
    sel = selection
    selection.on("contextmenu", function() {
      d3.event.preventDefault()
    })

    if (withSearch && d3.select("#mutation-search-input")[0][0] == null) {
      createMutationSearch(selection)
    }
    if (withToolbar && d3.select(".lineage-toolbar")[0][0] == null) {
      createToolBar(selection)
    }

    if (withTopScroll && d3.select(".top-scroll-wrapper")[0][0] == null) {
      selection
        .append("div")
        .attr("class", "top-scroll-wrapper")
        .append("div")
        .attr("class", "top-scroll-div")
    }

    if (
      selection[0][0].classList.contains("lineage") &&
      !d3.select(".lineage-wrapper")[0][0]
    ) {
      selection = selection.append("div").attr("class", "lineage-wrapper")
    }

    $(".top-scroll-wrapper").scroll(function() {
      $(".lineage-wrapper").scrollLeft($(".top-scroll-wrapper").scrollLeft())
    })
    $(".lineage-wrapper").scroll(function() {
      $(".top-scroll-wrapper").scrollLeft($(".lineage-wrapper").scrollLeft())
    })

    selection.each(function() {
      var minWidth = data.max_depth * minNodeWidth
      var containerWidth = d3
        .select(".lineage-wrapper")
        .node()
        .getBoundingClientRect().width
      var scrollWidth = Math.max(minWidth, containerWidth)
      data.x0 = height / 2
      data.y0 = 0
      if (vertical) {
        height = scrollWidth - margin.right - margin.left
        width = 80 + data.max_width * heightScalar
        nodeWidth = (1.4 * widthScalar * width) / data.max_depth
      } else {
        height = 80 + data.max_width * heightScalar
        width = scrollWidth - margin.right - margin.left
        nodeWidth = Math.max(
          minNodeWidth,
          (widthScalar * width) / data.max_depth
        )
      }
      tree.size([height - margin.top - margin.bottom, width])

      // Select the svg element, if it exists.
      var svg = selection.selectAll("svg").data([data])

      // Otherwise, create the skeletal chart.
      var svgEnter = svg
        .enter()
        .append("svg")
        .style("width", "100%")
      var svgDefs = svgEnter.append("defs")
      svgDefs
        .append("radialGradient")
        .attr("id", "unknown_gradient")
        .html(
          '<stop offset="10%" stop-color="#bcbcbc"/> <stop offset="80%" stop-color="#ccc"/> '
        )
      svgEnter.append("g");
      addDrawDropShadow(svg, dropShadow)

      var neededWidth = nodeWidth * data.max_depth + margin.right + margin.left
      d3.select(".top-scroll-div").style("min-width", neededWidth + "px")
      d3.select(".top-scroll-wrapper").style(
        "display",
        neededWidth <= containerWidth ? "none" : "block"
      )

      // Update the outer dimensions.
      svg.attr("height", height).style("min-width", neededWidth)

      // Update the inner dimensions.
      var g = svg
        .select("g")
        .attr("width", width)
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")")

      // assign parent, children, height, depth
      var root = data
      root.x0 = height / 2 // left edge of the rectangle
      root.y0 = 0 // top edge of the triangle

      update(root)

      function update(source) {
        // Compute the new tree layout.
        var nodes = tree.nodes(root).reverse(),
          links = tree.links(nodes)

        // Normalize for fixed-depth.
        nodes.forEach(function(d) {
          d.y = (d.depth - 1) * nodeWidth
        })

        if (vertical) {
          nodes.forEach(function(d) {
            var f = d.x
            d.x = d.y
            d.y = f
          })
        }

        // Update the nodes…
        var node = g.selectAll("g.node").data(nodes, function(d) {
          return d.id || (d.id = ++i)
        })

        // Enter any new nodes at the parent's previous position.
        var nodeEnter = node
          .enter()
          .append("g")
          .attr("class", function(d) {
            return (
              "node" +
              (d.hasOwnProperty("err") && d.err.length > 0 ? " has-err" : "")
            )
          })
          .attr("id", function(d) {
            return "node_" + d.slug
          })
          .attr("transform", function(d) {
            return "translate(" + source.y0 + "," + source.x0 + ")"
          })
          .on("contextmenu", click)

        nodeEnter
          .append("a")
          .attr("xlink:href", function(d) {
            return d.url
          })
          .append("circle")
          .attr("r", 1e-6)
          .style("filter", "url(#shadow);")
          .style("fill", function(d) {
            if (d.bg && d.bg.startsWith("linear:")) {
              svg
                .select("defs")
                .append("linearGradient")
                // the 0 -> XX is a hack to fix a weird name-changing bug
                .attr("id", d.slug.replace("0", "XX") + "_svggradient")
                .html(d.bg.replace("linear:", ""))
            }
          })
          .on("mouseover", function(d) {
            if (d.slug !== slug) {
              d3.select(this)
                .transition(150)
                .attr("r", function(d) {
                  return d._children ? defaultRadius : defaultRadius * 1.3
                })
            }

            const largeWindow = window.matchMedia("(min-width: 576px)").matches
            let dtext
            if (largeWindow) {
              dtext = `<strong>${d.name}</strong><br><span>`
            } else {
              dtext = `<strong><a href="${d.url}">${
                d.name
              }</a></strong><br><span>`
            }
            dtext += d.parent.name === "fakeroot" ? "" : d.parent.name
            if (d.mut) {
              let muts = d.mut.split("/")
              if (!show_inserts) {
                muts = muts.filter(d => (d.includes("ins") ? "" : d3))
                muts = muts.filter(d => (d.includes("ext") ? "" : d3))
              }
              if (!show_deletions) {
                muts = muts.filter(d => (d.includes("del") ? "" : d3))
              }
              dtext += ` &rarr; ${muts.join("/")}`
            }
            dtext += d.ref ? `<br><em>${d.ref}</em>` : ""
            dtext += "</span>"

            tooltip.html(dtext)

            if (largeWindow) {
              const _ttwidth = 200
              tooltip
                .style("width", _ttwidth + "px")
                .style("position", "absolute")
                .style("left", d3.event.pageX - _ttwidth / 2 + "px")
                .style("border-radius", "8px")
                .style("bottom", "inherit")
                .style("padding", ".6rem 0.5rem")
                .style("font-size", "inherit")
                .selectAll("span")
                .style("font-size", "0.75rem")
              tooltip.style(
                "top",
                d3.event.pageY - tooltip.node().clientHeight - 28 + "px"
              )
            } else {
              tooltip
                .style("width", "100%")
                .style("border-radius", "0px")
                .style("position", "fixed")
                .style("left", 0)
                .style("top", "inherit")
                .style("bottom", 0)
                .style("padding", "1rem 0.4rem 1.8rem 0.4rem")
                .style("font-size", "1.3rem")
                .selectAll("span")
                .style("font-size", "0.85rem")
              tooltip
                .append("div")
                .attr("class", "close-btn")
                .style("color", "#fff")
                .style("position", "absolute")
                .style("top", "10px")
                .style("right", "10px")
                .html("✖")
                .on("click", function() {
                  tooltip
                    .transition()
                    .duration(150)
                    .style("opacity", 0)
                    .transition()
                    .duration(0)
                    .style("left", -9999 + "px")
                })
            }
            tooltip
              .transition()
              .duration(150)
              .style("opacity", largeWindow ? 0.9 : 1)
          })
          .on("mouseout", function(d) {
            tooltip
              .transition()
              .duration(150)
              .style("opacity", 0)
              .transition()
              .duration(0)
              .style("left", -9999 + "px")

            if (d.slug !== slug) {
              d3.select(this)
                .transition(150)
                .attr("r", function(d) {
                  return d._children
                    ? defaultRadius / 2
                    : d.slug === slug
                    ? slugRadius
                    : defaultRadius
                })
            }
          })

        d3.select("#mutation-search-input").on("input", highlightMutations)
        d3.selectAll(".update-mutations").on("click", highlightMutations)

        //.style("fill", function(d) { return d._children ? "lightsteelblue" : "#fff"; });

        nodeEnter
          .append("text")
          .attr("x", function(d) {
            return text_position(d, slug, vertical)[0]
          })
          .attr("y", function(d) {
            return text_position(d, slug, vertical)[1]
          })
          .attr("dy", ".35em")
          .attr("text-anchor", function(d) {
            return text_position(d, slug, vertical)[2]
          })
          .text(function(d) {
            var t = d.name
            if (d.hasOwnProperty("err") && d.err.length > 0) {
              t += ` ! (${d.err[0].replace("SequenceMismatch:  diff: ", "")})`
            }
            return t
          })
          .attr("class", function(d) {
            return d.slug === slug ? "font-weight-bold" : ""
          })
          .style("fill-opacity", 1e-6)

        // Transition nodes to their new position.
        var nodeUpdate = node
          .transition()
          .duration(duration)
          .attr("transform", function(d) {
            return "translate(" + d.y + "," + d.x + ")"
          })

        nodeUpdate
          .select("circle")
          .attr("r", function(d) {
            return d._children
              ? defaultRadius / 2
              : d.slug === slug
              ? slugRadius
              : defaultRadius
          })
          .style("fill", function(d) {
            if (d.bg && d.bg.startsWith("linear:")) {
              return "url(#" + d.slug.replace("0", "XX") + "_svggradient)"
            } else if (d.bg === "?") {
              return "url(#unknown_gradient)"
            }
            return d.bg === "#222" ? "#888" : d.bg
          })
          .style("stroke-width", function(d) {
            return d._children ? defaultRadius / 2 + "px" : "1px"
          })

        nodeUpdate
          .select("text")
          .style("fill-opacity", 1)
          .attr("transform", function(d) {
            if (
              nodeWidth < 80 &&
              d.children &&
              d.children.length === 1 &&
              n_sibs(d) === 0 &&
              d.name.length > 8 &&
              n_sibs(d.parent) === 0
            ) {
              return "rotate(-15) translate(2, -2)"
            }
          })

        // Transition exiting nodes to the parent's new position.
        var nodeExit = node
          .exit()
          .transition()
          .duration(duration)
          .attr("transform", function(d) {
            return "translate(" + source.y + "," + source.x + ")"
          })
          .remove()

        nodeExit.select("circle").attr("r", 1e-6)

        nodeExit.select("text").style("fill-opacity", 1e-6)

        // Update the links…
        var link = g.selectAll("path.link").data(links, function(d) {
          return d.target.id
        })

        // Enter any new links at the parent's previous position.
        link
          .enter()
          .insert("path", "g")
          .attr("class", "link")
          .attr("d", function(d) {
            var o = { x: source.x0, y: source.y0 }
            return diagonal({ source: o, target: o })
          })

        // Transition links to their new position.
        link
          .transition()
          .duration(duration)
          .attr("d", diagonal)

        // Transition exiting nodes to the parent's new position.
        link
          .exit()
          .transition()
          .duration(duration)
          .attr("d", function(d) {
            var o = { x: source.x, y: source.y }
            return diagonal({ source: o, target: o })
          })
          .remove()

        //Option 1: remove node
        node.each(function(d) {
          if (d.name === "fakeroot") d3.select(this).remove()
        })

        link.each(function(d) {
          if (d.source.name === "fakeroot") d3.select(this).remove()
        })

        // Stash the old positions for transition.
        nodes.forEach(function(d) {
          d.x0 = d.x
          d.y0 = d.y
        })

        // toggle children on click
        function click(d) {
          if (d.children) {
            d._children = d.children
            d.children = null
          } else {
            d.children = d._children
            d._children = null
          }
          update(d)
        }

        function highlightMutations() {
          // these are hacks because the label class doesn't change to 'active'
          // fast enough
          let any
          let relparent
          if (this.classList.contains("mut-any")) {
            any = true
          } else if (this.classList.contains("mut-all")) {
            any = false
          } else {
            any = d3
              .select("#anytoggle")
              .node()
              .closest("label")
              .classList.contains("active")
          }
          if (this.classList.contains("mut-parent")) {
            relparent = true
          } else if (this.classList.contains("mut-root")) {
            relparent = false
          } else {
            relparent = d3
              .select("#parenttoggle")
              .node()
              .closest("label")
              .classList.contains("active")
          }

          var val = (d3.select("#mutation-search-input").node().value || "")
            .toUpperCase()
            .replace(",", " ")
            .split(" ")
            .filter(function(a) {
              return a.length > 1 ? a : null
            })

          if (val.length) {
            g.selectAll("path.link").attr("opacity", 0.35)
          } else {
            g.selectAll("path.link").attr("opacity", 1)
          }

          g.selectAll("circle")
            .attr("filter", function(d) {
              if (!val.length) {
                return null
              }
              if (any) {
                return val.some(function(v) {
                  return (relparent ? d.mut : d.rootmut).includes(v)
                })
                  ? "url(#dropshadow)"
                  : null
              }
              return val.every(function(v) {
                return (relparent ? d.mut : d.rootmut).includes(v)
              })
                ? "url(#dropshadow)"
                : null
            })
            .attr("opacity", function(d) {
              if (!val.length) {
                return 1
              }
              if (any) {
                return val.some(function(v) {
                  return (relparent ? d.mut : d.rootmut).includes(v)
                })
                  ? 1
                  : 0.3
              }
              return val.every(function(v) {
                return (relparent ? d.mut : d.rootmut).includes(v)
              })
                ? 1
                : 0.3
            })
          g.selectAll("text")
            .attr("opacity", function(d) {
              if (!val.length) {
                return 1
              }
              if (any) {
                return val.some(function(v) {
                  return (relparent ? d.mut : d.rootmut).includes(v)
                })
                  ? 1
                  : 0.3
              }
              return val.every(function(v) {
                return (relparent ? d.mut : d.rootmut).includes(v)
              })
                ? 1
                : 0.3
            })
            .style("font-weight", function(d) {
              if (!val.length) {
                return "inherit"
              }
              if (any) {
                return val.some(function(v) {
                  return (relparent ? d.mut : d.rootmut).includes(v)
                })
                  ? 500
                  : "inherit"
              }
              return val.every(function(v) {
                return (relparent ? d.mut : d.rootmut).includes(v)
              })
                ? 500
                : "inherit"
            })
        }
      }
    })
  }

  var resizeTimer
  $(window).on("resize", function(e) {
    clearTimeout(resizeTimer)
    resizeTimer = setTimeout(chart.update, 100)
  })

  // Public accessor methods

  chart.margin = function(_) {
    if (!arguments.length) return margin
    margin = _
    return chart
  }

  chart.width = function(value) {
    if (!arguments.length) return width
    width = value
    return chart
  }

  chart.slug = function(value) {
    if (!arguments.length) return slug
    slug = value
    return chart
  }

  chart.duration = function(value) {
    if (!arguments.length) return duration
    duration = value
    return chart
  }

  chart.height = function(value) {
    if (!arguments.length) return height
    height = value
    return chart
  }

  chart.heightScalar = function(value) {
    if (!arguments.length) return heightScalar
    heightScalar = value
    return chart
  }

  chart.scaleHeightUp = function() {
    heightScalar += 2
    chart(sel)
    return chart
  }

  chart.scaleHeightDown = function() {
    heightScalar -= 2
    chart(sel)
    return chart
  }

  chart.scaleWidthUp = function() {
    widthScalar += 0.07
    chart(sel)
    return chart
  }

  chart.scaleWidthDown = function() {
    widthScalar -= 0.07
    chart(sel)
    return chart
  }

  chart.widthScalar = function(value) {
    if (!arguments.length) return widthScalar
    widthScalar = value
    return chart
  }

  chart.withSearch = function(value) {
    if (!arguments.length) return withSearch
    withSearch = value
    return chart
  }

  chart.withToolbar = function(value) {
    if (!arguments.length) return withToolbar
    withToolbar = value
    return chart
  }

  function createToolBar(selection) {
    var tbar = selection
      .append("div")
      .attr({ class: "btn-toolbar lineage-toolbar", role: "toolbar" })
      .style("opacity", 0.8)
    var grp1 = tbar
      .append("div")
      .attr({ class: "btn-group btn-group-sm mr-2", role: "group" })
    grp1
      .append("button")
      .on("click", chart.scaleWidthDown)
      .attr({ type: "button", class: "btn btn-outline-dark" })
      .html("⇦")
    grp1
      .append("button")
      .on("click", chart.scaleWidthUp)
      .attr({ type: "button", class: "btn btn-outline-dark" })
      .html("⇨")
    grp1
      .append("button")
      .on("click", chart.scaleHeightDown)
      .attr({ type: "button", class: "btn btn-outline-dark" })
      .html("⇧")
    grp1
      .append("button")
      .on("click", chart.scaleHeightUp)
      .attr({ type: "button", class: "btn btn-outline-dark" })
      .html("⇩")
    var grp2 = tbar
      .append("div")
      .attr({ class: "btn-group btn-group-sm mr-2", role: "group" })
    grp2
      .append("button")
      .on("click", function() {
        chart.tree("tree")
      })
      .attr({ type: "button", class: "btn btn-outline-dark" })
      .html("⚯")
      .style("width", "2rem")
    grp2
      .append("button")
      .on("click", function() {
        chart.tree("cluster")
      })
      .attr({ type: "button", class: "btn btn-outline-dark" })
      .html("⚭")
      .style("width", "2rem")
  }

  chart.data = function(value) {
    if (!arguments.length) return data
    data = value
    // if (typeof updateData === "function") updateData();
    return chart
  }

  chart.update = function() {
    chart(sel)
  }

  chart.tree = function(value) {
    if (!arguments.length) return tree
    if (value === "cluster") {
      tree = d3.layout.cluster()
      chart(sel)
    } else {
      tree = d3.layout.tree()
      chart(sel)
    }
    return chart
  }

  return chart
}

/**
 * @see http://stackoverflow.com/questions/14865915/how-to-lower-the-opacity-of-the-alpha-layer-in-an-svg-filter/14871278#14871278
 */
function addDrawDropShadow(svg, dropShadow) {
  if (!d3.select("#dropshadow").node()) {
    var filter = svg
      .select("defs")
      .append("filter")
      .attr("id", "dropshadow")
      // x, y, width and height represent values in the current coordinate system that results
      // from taking the current user coordinate system in place at the time when the
      // <filter> element is referenced
      // (i.e., the user coordinate system for the element referencing the <filter> element via a filter attribute).
      .attr("filterUnits", "userSpaceOnUse")

    filter
      .append("feDropShadow")
      .attr("dx", dropShadow.dx || 0)
      .attr("dy", dropShadow.dy || 0)
      .attr("stdDeviation", dropShadow.stdDeviation || 4)
      .attr("flood-color", dropShadow.floodColor)
  }
}

function createMutationSearch(selection) {
  var wrapperDiv = selection
    .append("div")
    .append("div")
    .attr({ class: "row" })
  var searchDiv = wrapperDiv
    .append("div")
    .attr({ class: "input-group col-12 col-lg-8 mb-2" })
  searchDiv
    .append("div")
    .attr({ class: "input-group-prepend" })
    .append("span")
    .attr({ class: "input-group-text" })
    .text("Search")

  searchDiv.append("input").attr({
    type: "search",
    class: "form-control",
    name: "textInput",
    placeholder: "Mutations (e.g. A206K) separated by spaces",
    id: "mutation-search-input"
  })

  var btngroup = searchDiv.append("div").attr({ class: "input-group-append" })

  var anyallgroup = btngroup
    .append("div")
    .attr({ class: "btn-group-toggle btn-group", "data-toggle": "buttons" })

  anyallgroup
    .append("label")
    .attr({ class: "btn btn-outline-primary update-mutations mut-all active" })
    .text("all")
    .append("input")
    .attr({
      type: "radio",
      name: "anyall",
      id: "alltoggle",
      autocomplete: "off"
    })

  anyallgroup
    .append("label")
    .attr({ class: "btn btn-outline-primary update-mutations mut-any" })
    .text("any")
    .append("input")
    .attr({
      type: "radio",
      name: "anyall",
      id: "anytoggle",
      autocomplete: "off"
    })

  var rightdiv = wrapperDiv
    .append("div")
    .attr({ class: "input-group col-12 col-lg-4 mb-2" })
  rightdiv
    .append("div")
    .attr({ class: "input-group-prepend" })
    .append("span")
    .attr({ class: "input-group-text" })
    .text("Relative to")

  var relativetogroup = rightdiv.append("div").attr({
    class: "btn-group-toggle btn-group input-group-append",
    "data-toggle": "buttons"
  })

  relativetogroup
    .append("label")
    .attr({ class: "btn btn-outline-primary update-mutations mut-root active" })
    .text("root")
    .append("input")
    .attr({
      type: "radio",
      name: "parentroot",
      id: "roottoggle",
      autocomplete: "off"
    })

  relativetogroup
    .append("label")
    .attr({ class: "btn btn-outline-primary update-mutations mut-parent" })
    .text("parent")
    .append("input")
    .attr({
      type: "radio",
      name: "parentroot",
      id: "parenttoggle",
      autocomplete: "off"
    })
}

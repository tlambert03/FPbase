// Use modular D3 imports for tree-shaking (only in d3Charts bundle)

import { max, min } from "d3-array"
import { cluster, hierarchy, tree } from "d3-hierarchy"
import { select, selectAll } from "d3-selection"

const $ = window.jQuery // jQuery loaded from CDN

// Create d3 namespace for backward compatibility with existing code
const d3 = { select, selectAll, hierarchy, tree, cluster, max, min }

function n_sibs(node) {
  return node.parent ? node.parent.children.length - 1 : 0
}

function is_first_sib(node) {
  if (!node.parent || node.parent.children.length === 1) {
    return false
  }
  const i = node.parent.children.map((e) => e.name).indexOf(node.name)
  if (i === 0) {
    return true
  }
  return false
}

function is_last_sib(node) {
  if (!node.parent || node.parent.children.length === 1) {
    return false
  }
  const i = node.parent.children.map((e) => e.name).indexOf(node.name)
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
  const config = conf || {}
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
    tree = d3.tree(),
    withTopScroll = config.withTopScroll || false,
    withSearch = config.withSearch || false,
    withToolbar = config.withToolbar || false,
    vertical = config.vertical || false,
    dropShadow = {
      stdDeviation: 4,
      dx: 0,
      dy: 0,
      floodColor: "#38e",
    }

  // Custom diagonal function for D3 v7 (replaces d3.svg.diagonal)
  var diagonal = (d) => `M${d.source.y},${d.source.x}
            C${(d.source.y + d.target.y) / 2},${d.source.x}
             ${(d.source.y + d.target.y) / 2},${d.target.x}
             ${d.target.y},${d.target.x}`

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
    selection.on("contextmenu", (event) => {
      event.preventDefault()
    })

    if (withSearch && d3.select("#mutation-search-input").node() == null) {
      createMutationSearch(selection)
    }
    if (withToolbar && d3.select(".lineage-toolbar").node() == null) {
      createToolBar(selection)
    }

    if (withTopScroll && d3.select(".top-scroll-wrapper").node() == null) {
      selection
        .append("div")
        .attr("class", "top-scroll-wrapper")
        .append("div")
        .attr("class", "top-scroll-div")
    }

    if (selection.node().classList.contains("lineage") && !d3.select(".lineage-wrapper").node()) {
      selection = selection.append("div").attr("class", "lineage-wrapper")
    }

    $(".top-scroll-wrapper").scroll(() => {
      $(".lineage-wrapper").scrollLeft($(".top-scroll-wrapper").scrollLeft())
    })
    $(".lineage-wrapper").scroll(() => {
      $(".top-scroll-wrapper").scrollLeft($(".lineage-wrapper").scrollLeft())
    })

    selection.each(() => {
      var minWidth = data.max_depth * minNodeWidth
      var containerWidth = d3.select(".lineage-wrapper").node().getBoundingClientRect().width
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
        nodeWidth = Math.max(minNodeWidth, (widthScalar * width) / data.max_depth)
      }
      tree.size([height - margin.top - margin.bottom, width])

      // Select the svg element, if it exists.
      var svg = selection.selectAll("svg").data([data])

      // Otherwise, create the skeletal chart.
      var svgEnter = svg.enter().append("svg").style("width", "100%")
      var svgDefs = svgEnter.append("defs")
      svgDefs
        .append("radialGradient")
        .attr("id", "unknown_gradient")
        .html('<stop offset="10%" stop-color="#bcbcbc"/> <stop offset="80%" stop-color="#ccc"/> ')
      svgEnter.append("g")
      addDrawDropShadow(svg, dropShadow)

      var neededWidth = nodeWidth * data.max_depth + margin.right + margin.left
      d3.select(".top-scroll-div").style("min-width", `${neededWidth}px`)
      d3.select(".top-scroll-wrapper").style(
        "display",
        neededWidth <= containerWidth ? "none" : "block"
      )

      // Merge enter and update selections for svg
      svg = svg.merge(svgEnter)

      // Update the outer dimensions on the merged selection
      svg.attr("height", height).style("min-width", neededWidth)

      // Update the inner dimensions.
      var g = svg
        .select("g")
        .attr("width", width)
        .attr("transform", `translate(${margin.left},${margin.top})`)

      // D3 v7: Create hierarchy from data if not already a hierarchy
      // Check if data is already a d3 hierarchy (has descendants method)
      var root = data.descendants ? data : d3.hierarchy(data)
      root.x0 = height / 2 // left edge of the rectangle
      root.y0 = 0 // top edge of the triangle

      update(root)

      function update(source) {
        // Compute the new tree layout.
        // D3 v7: tree() returns the hierarchy, then we get nodes and links
        var hierarchy = tree(root)
        var nodes = hierarchy.descendants().reverse()
        var links = hierarchy.links()

        // Normalize for fixed-depth.
        nodes.forEach((d) => {
          d.y = (d.depth - 1) * nodeWidth
        })

        // Check for negative Y positions and shift if needed
        // But exclude fakeroot from this calculation since it will be removed from display
        var visibleNodes = nodes.filter((d) => d.data.name !== "fakeroot")
        const minY = d3.min(visibleNodes, (d) => d.y)
        if (minY < 0) {
          const yShift = -minY
          nodes.forEach((d) => {
            d.y = d.y + yShift
          })
        }

        if (vertical) {
          nodes.forEach((d) => {
            var f = d.x
            d.x = d.y
            d.y = f
          })
        }

        // Calculate actual height needed based on node positions
        var minX = d3.min(nodes, (d) => d.x)
        var maxX = d3.max(nodes, (d) => d.x)
        var neededHeight = maxX - minX + margin.top + margin.bottom + 40 // extra padding

        // Update SVG height if needed
        if (neededHeight > height) {
          height = neededHeight
          svg.attr("height", height)
        }

        // Update the nodes…
        var node = g.selectAll("g.node").data(nodes, (d) => {
          if (!d.id) {
            d.id = ++i
          }
          return d.id
        })

        // Enter any new nodes at the parent's previous position.
        var nodeEnter = node
          .enter()
          .append("g")
          .attr("class", (d) => `node${d.data.err && d.data.err.length > 0 ? " has-err" : ""}`)
          .attr("id", (d) => `node_${d.data.slug}`)
          .attr("transform", (_d) => `translate(${source.y0},${source.x0})`)
          .on("contextmenu", click)

        nodeEnter
          .append("a")
          .attr("xlink:href", (d) => d.data.url)
          .append("circle")
          .attr("r", 1e-6)
          .style("filter", "url(#shadow);")
          .style("fill", (d) => {
            if (d.data.bg?.startsWith("linear:")) {
              svg
                .select("defs")
                .append("linearGradient")
                // the 0 -> XX is a hack to fix a weird name-changing bug
                .attr("id", `${d.data.slug.replace("0", "XX")}_svggradient`)
                .html(d.data.bg.replace("linear:", ""))
            }
          })
          .on("mouseover", function (event, d) {
            if (d.data.slug !== slug) {
              d3.select(this)
                .transition(150)
                .attr("r", (d) => (d._children ? defaultRadius : defaultRadius * 1.3))
            }

            const largeWindow = window.matchMedia("(min-width: 576px)").matches
            let dtext
            if (largeWindow) {
              dtext = `<strong>${d.data.name}</strong><br><span>`
            } else {
              dtext = `<strong><a href="${d.data.url}">${d.data.name}</a></strong><br><span>`
            }
            dtext += d.parent ? (d.parent.data.name === "fakeroot" ? "" : d.parent.data.name) : ""
            if (d.data.mut) {
              let muts = d.data.mut.split("/")
              if (!show_inserts) {
                muts = muts.filter((m) => (m.includes("ins") ? "" : m))
                muts = muts.filter((m) => (m.includes("ext") ? "" : m))
              }
              if (!show_deletions) {
                muts = muts.filter((m) => (m.includes("del") ? "" : m))
              }
              dtext += ` &rarr; ${muts.join("/")}`
            }
            dtext += d.data.ref ? `<br><em>${d.data.ref}</em>` : ""
            dtext += "</span>"

            tooltip.html(dtext)

            if (largeWindow) {
              const _ttwidth = 200
              tooltip
                .style("width", `${_ttwidth}px`)
                .style("position", "absolute")
                .style("left", `${event.pageX - _ttwidth / 2}px`)
                .style("border-radius", "8px")
                .style("bottom", "inherit")
                .style("padding", ".6rem 0.5rem")
                .style("font-size", "inherit")
                .selectAll("span")
                .style("font-size", "0.75rem")
              tooltip.style("top", `${event.pageY - tooltip.node().clientHeight - 28}px`)
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
                .on("click", () => {
                  tooltip
                    .transition()
                    .duration(150)
                    .style("opacity", 0)
                    .transition()
                    .duration(0)
                    .style("left", `${-9999}px`)
                })
            }
            tooltip
              .transition()
              .duration(150)
              .style("opacity", largeWindow ? 0.9 : 1)
          })
          .on("mouseout", function (_event, d) {
            tooltip
              .transition()
              .duration(150)
              .style("opacity", 0)
              .transition()
              .duration(0)
              .style("left", `${-9999}px`)

            if (d.data.slug !== slug) {
              d3.select(this)
                .transition(150)
                .attr("r", (d) =>
                  d._children
                    ? defaultRadius / 2
                    : d.data.slug === slug
                      ? slugRadius
                      : defaultRadius
                )
            }
          })

        d3.select("#mutation-search-input").on("input", highlightMutations)
        d3.selectAll(".update-mutations").on("click", highlightMutations)

        //.style("fill", function(d) { return d._children ? "lightsteelblue" : "#fff"; });

        nodeEnter
          .append("text")
          .attr("x", (d) => text_position(d, slug, vertical)[0])
          .attr("y", (d) => text_position(d, slug, vertical)[1])
          .attr("dy", ".35em")
          .attr("text-anchor", (d) => text_position(d, slug, vertical)[2])
          .text((d) => {
            var t = d.data.name
            if (d.data.err && d.data.err.length > 0) {
              t += ` ! (${d.data.err[0].replace("SequenceMismatch:  diff: ", "")})`
            }
            return t
          })
          .attr("class", (d) => (d.data.slug === slug ? "font-weight-bold" : ""))
          .style("fill-opacity", 1e-6)

        // Merge enter and update selections
        var nodeUpdate = node
          .merge(nodeEnter)
          .transition()
          .duration(duration)
          .attr("transform", (d) => `translate(${d.y},${d.x})`)

        nodeUpdate
          .select("circle")
          .attr("r", (d) =>
            d._children ? defaultRadius / 2 : d.data.slug === slug ? slugRadius : defaultRadius
          )
          .style("fill", (d) => {
            if (d.data.bg?.startsWith("linear:")) {
              return `url(#${d.data.slug.replace("0", "XX")}_svggradient)`
            } else if (d.data.bg === "?") {
              return "url(#unknown_gradient)"
            }
            return d.data.bg === "#222" ? "#888" : d.data.bg
          })
          .style("stroke-width", (d) => (d._children ? `${defaultRadius / 2}px` : "1px"))

        nodeUpdate
          .select("text")
          .style("fill-opacity", 1)
          .attr("transform", (d) => {
            if (
              nodeWidth < 80 &&
              d.children &&
              d.children.length === 1 &&
              n_sibs(d) === 0 &&
              d.data.name.length > 8 &&
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
          .attr("transform", (_d) => `translate(${source.y},${source.x})`)
          .remove()

        nodeExit.select("circle").attr("r", 1e-6)

        nodeExit.select("text").style("fill-opacity", 1e-6)

        // Update the links…
        var link = g.selectAll("path.link").data(links, (d) => d.target.id)

        // Enter any new links at the parent's previous position.
        var linkEnter = link
          .enter()
          .insert("path", "g")
          .attr("class", "link")
          .attr("d", (_d) => {
            var o = { x: source.x0, y: source.y0 }
            return diagonal({ source: o, target: o })
          })

        // Transition links to their new position.
        link.merge(linkEnter).transition().duration(duration).attr("d", diagonal)

        // Transition exiting nodes to the parent's new position.
        link
          .exit()
          .transition()
          .duration(duration)
          .attr("d", (_d) => {
            var o = { x: source.x, y: source.y }
            return diagonal({ source: o, target: o })
          })
          .remove()

        // Remove fakeroot nodes and links
        // Must operate on merged selection to catch both enter and update selections
        node.merge(nodeEnter).each(function (d) {
          if (d.data.name === "fakeroot") d3.select(this).remove()
        })

        link.merge(linkEnter).each(function (d) {
          if (d.source.data.name === "fakeroot") d3.select(this).remove()
        })

        // Stash the old positions for transition.
        nodes.forEach((d) => {
          d.x0 = d.x
          d.y0 = d.y
        })

        // toggle children on click
        function click(_event, d) {
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
            any = d3.select("#anytoggle").node().closest("label").classList.contains("active")
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
            .filter((a) => (a.length > 1 ? a : null))

          if (val.length) {
            g.selectAll("path.link").attr("opacity", 0.35)
          } else {
            g.selectAll("path.link").attr("opacity", 1)
          }

          g.selectAll("circle")
            .attr("filter", (d) => {
              if (!val.length) {
                return null
              }
              if (any) {
                return val.some((v) => (relparent ? d.data.mut : d.data.rootmut).includes(v))
                  ? "url(#dropshadow)"
                  : null
              }
              return val.every((v) => (relparent ? d.data.mut : d.data.rootmut).includes(v))
                ? "url(#dropshadow)"
                : null
            })
            .attr("opacity", (d) => {
              if (!val.length) {
                return 1
              }
              if (any) {
                return val.some((v) => (relparent ? d.data.mut : d.data.rootmut).includes(v))
                  ? 1
                  : 0.3
              }
              return val.every((v) => (relparent ? d.data.mut : d.data.rootmut).includes(v))
                ? 1
                : 0.3
            })
          g.selectAll("text")
            .attr("opacity", (d) => {
              if (!val.length) {
                return 1
              }
              if (any) {
                return val.some((v) => (relparent ? d.data.mut : d.data.rootmut).includes(v))
                  ? 1
                  : 0.3
              }
              return val.every((v) => (relparent ? d.data.mut : d.data.rootmut).includes(v))
                ? 1
                : 0.3
            })
            .style("font-weight", (d) => {
              if (!val.length) {
                return "inherit"
              }
              if (any) {
                return val.some((v) => (relparent ? d.data.mut : d.data.rootmut).includes(v))
                  ? 500
                  : "inherit"
              }
              return val.every((v) => (relparent ? d.data.mut : d.data.rootmut).includes(v))
                ? 500
                : "inherit"
            })
        }
      }
    })
  }

  var resizeTimer
  $(window).on("resize", (_e) => {
    clearTimeout(resizeTimer)
    resizeTimer = setTimeout(chart.update, 100)
  })

  // Public accessor methods

  chart.margin = (...args) => {
    if (!args.length) return margin
    margin = args[0]
    return chart
  }

  chart.width = (...args) => {
    if (!args.length) return width
    width = args[0]
    return chart
  }

  chart.slug = (...args) => {
    if (!args.length) return slug
    slug = args[0]
    return chart
  }

  chart.duration = (...args) => {
    if (!args.length) return duration
    duration = args[0]
    return chart
  }

  chart.height = (...args) => {
    if (!args.length) return height
    height = args[0]
    return chart
  }

  chart.heightScalar = (...args) => {
    if (!args.length) return heightScalar
    heightScalar = args[0]
    return chart
  }

  chart.scaleHeightUp = () => {
    heightScalar += 2
    chart(sel)
    return chart
  }

  chart.scaleHeightDown = () => {
    heightScalar -= 2
    chart(sel)
    return chart
  }

  chart.scaleWidthUp = () => {
    widthScalar += 0.07
    chart(sel)
    return chart
  }

  chart.scaleWidthDown = () => {
    widthScalar -= 0.07
    chart(sel)
    return chart
  }

  chart.widthScalar = (...args) => {
    if (!args.length) return widthScalar
    widthScalar = args[0]
    return chart
  }

  chart.withSearch = (...args) => {
    if (!args.length) return withSearch
    withSearch = args[0]
    return chart
  }

  chart.withToolbar = (...args) => {
    if (!args.length) return withToolbar
    withToolbar = args[0]
    return chart
  }

  function createToolBar(selection) {
    var tbar = selection
      .append("div")
      .attr("class", "btn-toolbar lineage-toolbar")
      .attr("role", "toolbar")
      .style("opacity", 0.8)
    var grp1 = tbar.append("div").attr("class", "btn-group btn-group-sm mr-2").attr("role", "group")
    grp1
      .append("button")
      .on("click", chart.scaleWidthDown)
      .attr("type", "button")
      .attr("class", "btn btn-outline-dark")
      .html("⇦")
    grp1
      .append("button")
      .on("click", chart.scaleWidthUp)
      .attr("type", "button")
      .attr("class", "btn btn-outline-dark")
      .html("⇨")
    grp1
      .append("button")
      .on("click", chart.scaleHeightDown)
      .attr("type", "button")
      .attr("class", "btn btn-outline-dark")
      .html("⇧")
    grp1
      .append("button")
      .on("click", chart.scaleHeightUp)
      .attr("type", "button")
      .attr("class", "btn btn-outline-dark")
      .html("⇩")
    var grp2 = tbar.append("div").attr("class", "btn-group btn-group-sm mr-2").attr("role", "group")
    grp2
      .append("button")
      .on("click", () => {
        chart.tree("tree")
      })
      .attr("type", "button")
      .attr("class", "btn btn-outline-dark")
      .html("⚯")
      .style("width", "2rem")
    grp2
      .append("button")
      .on("click", () => {
        chart.tree("cluster")
      })
      .attr("type", "button")
      .attr("class", "btn btn-outline-dark")
      .html("⚭")
      .style("width", "2rem")
  }

  chart.data = (...args) => {
    if (!args.length) return data
    data = args[0]
    // if (typeof updateData === "function") updateData();
    return chart
  }

  chart.update = () => {
    chart(sel)
  }

  chart.tree = (...args) => {
    if (!args.length) return tree
    if (args[0] === "cluster") {
      tree = d3.cluster()
      chart(sel)
    } else {
      tree = d3.tree()
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
    const filter = svg
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
  var wrapperDiv = selection.append("div").append("div").attr("class", "row")
  var searchDiv = wrapperDiv.append("div").attr("class", "input-group col-12 col-lg-8 mb-2")
  searchDiv
    .append("div")
    .attr("class", "input-group-prepend")
    .append("span")
    .attr("class", "input-group-text")
    .text("Search")

  searchDiv
    .append("input")
    .attr("type", "search")
    .attr("class", "form-control")
    .attr("name", "textInput")
    .attr("placeholder", "Mutations (e.g. A206K) separated by spaces")
    .attr("id", "mutation-search-input")

  var btngroup = searchDiv.append("div").attr("class", "input-group-append")

  var anyallgroup = btngroup
    .append("div")
    .attr("class", "btn-group-toggle btn-group")
    .attr("data-toggle", "buttons")

  anyallgroup
    .append("label")
    .attr("class", "btn btn-outline-primary update-mutations mut-all active")
    .text("all")
    .append("input")
    .attr("type", "radio")
    .attr("name", "anyall")
    .attr("id", "alltoggle")
    .attr("autocomplete", "off")

  anyallgroup
    .append("label")
    .attr("class", "btn btn-outline-primary update-mutations mut-any")
    .text("any")
    .append("input")
    .attr("type", "radio")
    .attr("name", "anyall")
    .attr("id", "anytoggle")
    .attr("autocomplete", "off")

  var rightdiv = wrapperDiv.append("div").attr("class", "input-group col-12 col-lg-4 mb-2")
  rightdiv
    .append("div")
    .attr("class", "input-group-prepend")
    .append("span")
    .attr("class", "input-group-text")
    .text("Relative to")

  var relativetogroup = rightdiv
    .append("div")
    .attr("class", "btn-group-toggle btn-group input-group-append")
    .attr("data-toggle", "buttons")

  relativetogroup
    .append("label")
    .attr("class", "btn btn-outline-primary update-mutations mut-root active")
    .text("root")
    .append("input")
    .attr("type", "radio")
    .attr("name", "parentroot")
    .attr("id", "roottoggle")
    .attr("autocomplete", "off")

  relativetogroup
    .append("label")
    .attr("class", "btn btn-outline-primary update-mutations mut-parent")
    .text("parent")
    .append("input")
    .attr("type", "radio")
    .attr("name", "parentroot")
    .attr("id", "parenttoggle")
    .attr("autocomplete", "off")
}
// Force rebuild

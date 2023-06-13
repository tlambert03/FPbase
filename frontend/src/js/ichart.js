import noUiSlider from "nouislider"
import d3 from "d3"
import $ from "jquery"

export default function FPPropChart() {
  const margin = { top: 20, right: 30, bottom: 20, left: 68 }
  const chartWidth = 700
  const chartHeight = 700
  let width = chartWidth - margin.right
  let height = chartHeight - margin.top - margin.bottom
  // sel,
  let currentData
  let FPdata = [] // Where the fluorescent protein data table will end up.;
  let currentX = "ex_max"
  let currentY = "em_max"
  // danceInterval,
  const circleMoveDuration = 400
  const symbolsize = 8 // radius of circle
  const bigscale = 1.5 // how much to scale up on mouseover

  width = chartWidth - margin.right
  height = chartHeight - margin.top - margin.bottom

  // global variable to set the ranges over which the data is filtered.
  const dataFilters = {
    ex_max: [350, 800, 1], // array values represent [min range, max range, step (for the range slider)]
    em_max: [350, 800, 1],
    ext_coeff: [10000, 230000, 1000],
    qy: [0, 1, 0.01],
    brightness: [0, 165, 1],
    agg: "",
  }
  // string variables for updating the axis labels
  const strings = {
    em_max: "Emission Wavelength (nm)",
    ex_max: "Excitation Wavelength (nm)",
    stokes: "Stokes Shift (nm)",
    ext_coeff: "Extinction Coefficient",
    qy: "Quantum Yield",
    brightness: "Brightness",
    pka: "pKa",
    bleach: "Bleaching Half-life (s)",
    maturation: "Maturation Half-time (min)",
    lifetime: "Lifetime (ns)",
  }

  // Scales and axes
  const xScale = d3.scale.linear().range([0, width])

  const yScale = d3.scale.linear().range([height, 0])

  // This scale will set the saturation (gray to saturated color).  We will use it for mapping brightness.
  const saturationScale = d3.scale
    .linear()
    .range([0, 1])
    .domain([0, 100])

  // This scale will set the hue.  We will use it for mapping emission wavelength.
  const hueScale = d3.scale
    .linear()
    .range([300, 300, 240, 0, 0])
    .domain([200, 405, 440, 650, 850])

  // X and Y axes
  const xAxisBottom = d3.svg
    .axis()
    .scale(xScale)
    .tickSize(5)
    .tickSubdivide(true)
  const yAxisLeft = d3.svg
    .axis()
    .scale(yScale)
    .tickSize(5)
    .orient("left")
    .tickSubdivide(true)

  // top and right axes are identical but without tick labels
  const xAxisTop = d3.svg
    .axis()
    .scale(xScale)
    .tickSize(5)
    .orient("top")
    .tickSubdivide(true)
    .tickFormat(function(d) {
      return ""
    })
  const yAxisRight = d3.svg
    .axis()
    .scale(yScale)
    .tickSize(5)
    .orient("right")
    .tickSubdivide(true)
    .tickFormat(function(d) {
      return ""
    })

  // on page load, listen to slider events and respond by updating the filter ranges (and updating the ui)

  function chart(selection) {
    // sel = selection

    // Create the SVG container and set the origin.

    // Select the svg element, if it exists.
    let svg = selection.selectAll("svg").data([FPdata])

    // Otherwise, create the skeletal chart.
    const svgEnter = svg.enter().append("svg")
    const gEnter = svgEnter.append("g")

    function drawGraph() {
      // redraw axes with new domainsc
      svg.select(".x.axis.bottom").call(xAxisBottom)
      svg.select(".y.axis.left").call(yAxisLeft)
      svg.select(".x.axis.top").call(xAxisTop)
      svg.select(".y.axis.right").call(yAxisRight)

      svg
        .selectAll("circle.FP")
        .attr("cx", function(d) {
          return xScale(d[currentX])
        })
        .attr("cy", function(d) {
          return yScale(d[currentY])
        })

      svg
        .selectAll("rect.FP")
        .attr("x", function(d) {
          return xScale(d[currentX]) - symbolsize
        })
        .attr("y", function(d) {
          return yScale(d[currentY]) - symbolsize
        })

      svg
        .selectAll("text.FP")
        .attr("x", function(d) {
          return xScale(d[currentX]) - symbolsize / 2
        })
        .attr("y", function(d) {
          return yScale(d[currentY]) + symbolsize / 2
        })
    }

    svg
      .attr("viewBox", "0 0 760 760")
      .attr("id", "mainchart")
      .attr("preserveAspectRatio", "xMinYMin meet")

    svg = svg.select("g")

    gEnter
      .attr("transform", `translate(${margin.left},${margin.top})`)
      .classed("svg-content-responsive", true)

    // Add the axes
    gEnter
      .append("g")
      .attr("class", "x axis bottom")
      .attr("transform", `translate(0,${height})`)
      .call(xAxisBottom)
    gEnter
      .append("g")
      .attr("class", "y axis left")
      .call(yAxisLeft)
    gEnter
      .append("g")
      .attr("class", "x axis top")
      .call(xAxisTop)
    gEnter
      .append("svg:g")
      .attr("class", "y axis right")
      .attr("transform", `translate(${width},0)`)
      .call(yAxisRight)

    // Add an x-axis label.
    gEnter
      .append("text")
      .attr("class", "x label")
      .attr("text-anchor", "middle")
      .attr("x", width / 2)
      .attr("y", height - 10)
      .text("Excitation wavelength (nm)")

    // Add a y-axis label.
    gEnter
      .append("text")
      .attr("class", "y label")
      .attr("text-anchor", "middle")
      .attr("x", -height / 2)
      .attr("y", margin.left - 30)
      .attr("transform", "rotate(-90)")
      .text("Emission wavelength (nm)")

    // Add a clipping path so that data points don't go outside of frame
    gEnter
      .append("clipPath") // Make a new clipPath
      .attr("id", "chart-area") // Assign an ID
      .append("rect")
      .attr("width", width)
      .attr("height", height)

    // enable zooming
    const zoom = d3.behavior
      .zoom()
      .x(xScale)
      .y(yScale)
      .scaleExtent([1, 10])
      .on("zoom", drawGraph)

    gEnter
      .append("rect")
      .attr("class", "pane")
      .attr("width", width)
      .attr("height", height)
      .call(zoom)

    function addactions(sel) {
      sel.on("click", function(e) {
        window.location = e.url
      })
      sel
        .on("mouseover", function(d) {
          // Get this bar's x/y values, then augment for the tooltip
          if (d3.select(this).attr("cx")) {
            // if circle
            d3.select(this)
              .transition()
              .duration(100)
              .attr("r", symbolsize * bigscale)
          } else if (d3.select(this).attr("x")) {
            // if rectangle
            d3.select(this)
              .transition()
              .duration(100)
              .attr("x", function(d_) {
                return xScale(d_[currentX]) - symbolsize * bigscale
              })
              .attr("y", function(d_) {
                return yScale(d_[currentY]) - symbolsize * bigscale
              })
              .attr("width", symbolsize * 2 * bigscale)
              .attr("height", symbolsize * 2 * bigscale)
          }
          const ypos = d3.event.pageY - $("rect.pane").position().top
          const xpos = d3.event.pageX - $("rect.pane").position().left

          const yPosition = $("rect.pane").position().top + 10
          let xPosition
          if (xpos < 179 && ypos < 141) {
            xPosition =
              $("#mainchart").position().left + $("#mainchart").width() - 180
          } else {
            xPosition = $("rect.pane").position().left + 10
          }
          // Update the tooltip position and value
          d3.select("#tooltip")
            .style("left", `${xPosition}px`)
            .style("top", `${yPosition}px`)
            .select("#exvalue")
            .text(d.ex_max)
          d3.select("#tooltip")
            .select("#emvalue")
            .text(d.em_max)
          d3.select("#tooltip")
            .select("#ecvalue")
            .text(d.ext_coeff)
          d3.select("#tooltip")
            .select("#qyvalue")
            .text(d.qy)
          d3.select("#tooltip")
            .select("h3")
            .html(d.name)
          d3.select("#tooltip")
            .select("#brightnessvalue")
            .text(d.brightness)

          // Show the tooltip
          d3.select("#tooltip").classed("hidden", false)
        })

        .on("mouseout", function() {
          if (d3.select(this).attr("cx")) {
            // if circle
            d3.select(this)
              .transition()
              .duration(200)
              .attr("r", symbolsize)
          } else if (d3.select(this).attr("x")) {
            // if circle
            d3.select(this)
              .transition()
              .duration(200)
              .attr("x", function(d) {
                return xScale(d[currentX]) - symbolsize
              })
              .attr("y", function(d) {
                return yScale(d[currentY]) - symbolsize
              })
              .attr("width", symbolsize * 2)
              .attr("height", symbolsize * 2)
          }
          // Hide the tooltip
          d3.select("#tooltip").classed("hidden", true)
        })
    }

    function plotcircle(sel) {
      const circle = sel
        .append("circle")
        .attr("class", "FP")
        .attr("r", symbolsize)
        .attr("stroke", "#000")
        .attr("opacity", 0.7)
        .style("fill", function(d) {
          return d3.hsl(hueScale(d.em_max), saturationScale(d.brightness), 0.5)
        })
      addactions(circle)
    }

    function plotsquare(sel) {
      const square = sel
        .append("rect")
        .attr("class", "FP")
        .attr("width", symbolsize * 2)
        .attr("height", symbolsize * 2)
        .attr("stroke", "#000")
        .attr("opacity", 0.7)
        .style("fill", function(d) {
          return d3.hsl(hueScale(d.em_max), saturationScale(d.brightness), 0.5)
        })
      addactions(square)
    }

    function labelAgg(sel) {
      sel
        .append("text")
        .attr("class", "FP")
        .text(function(d) {
          if (d.agg === "d") {
            return "2"
          }
          if (d.agg === "td") {
            return "t"
          }
          if (d.agg === "wd") {
            return "2"
          }
          if (d.agg === "t") {
            return "4"
          }
          return ""
        })
    }

    // i added this more flexible plotting function to be able to plot different variables on each axis.
    // It takes three optional parameters: the data array, and two axes variables.
    function plot() {
      // set default values... if plot() is called without arguments, these default values will be used.
      const xvar = currentX
      const yvar = currentY
      currentData = FPdata

      // helper function to iterate through all of the data filters (without having to type them all out)
      function filtercheck(datum) {
        // eslint-disable-next-line
        for (const f in dataFilters) {
          if (f === "agg") {
            if (dataFilters[f]) {
              if (datum[f] !== dataFilters[f]) {
                return false
              }
            }
          } else {
            const v = dataFilters[f]
            if (datum[f] < v[0] || datum[f] > v[1]) {
              return false
            }
          }
        }
        return true
      }

      // filter the currentData according to the user settings for EC, QY, and brightness range
      currentData = currentData.filter(function(d) {
        return filtercheck(d) ? d : null
      })

      // filter out currentData with empty values
      currentData = currentData.filter(function(d) {
        return d[xvar] > 0 && d[yvar] > 0
      })

      // update scale domains based on currentData
      xScale
        .domain([
          d3.min(currentData, function(d) {
            return 0.99 * d[xvar]
          }),
          d3.max(currentData, function(d) {
            return 1.01 * d[xvar]
          }),
        ])
        .nice()
      zoom.x(xScale)

      yScale
        .domain([
          d3.min(currentData, function(d) {
            return 0.99 * d[yvar]
          }),
          d3.max(currentData, function(d) {
            return 1.01 * d[yvar]
          }),
        ])
        .nice()
      zoom.y(yScale)

      // relabel X and Y axes
      svg.select(".x.label").text(strings[xvar])
      svg.select(".y.label").text(strings[yvar])

      // Join new currentData with old elements, if any.
      const datagroup = svg.selectAll("g.FP").data(currentData, function(d) {
        return d.name
      })
      const entergroup = datagroup
        .enter()
        .append("g")
        .attr("class", "FP")
        .attr("clip-path", "url(#chart-area)")
        .call(zoom) // so we can zoom while moused over elements

      entergroup.each(function(d, i) {
        // determine type of protein and whether to plot a circle or a square
        if (d.cofactor !== "") {
          // plot squeares (for proteins with cofactor)
          plotsquare(d3.select(this))
        } else {
          // plot new squares
          plotcircle(d3.select(this))
        }
        // add text to markers
        labelAgg(d3.select(this))
      })

      // Remove old elements as needed.
      datagroup.exit().remove()

      // move circles to their new positions (based on axes) with transition animation
      datagroup.each(function(d, i) {
        const current = d3.select(this)
        current
          .selectAll("circle.FP")
          .transition()
          .attr("cx", function(d_) {
            return xScale(d_[xvar])
          })
          .attr("cy", function(d_) {
            return yScale(d_[yvar])
          })
          .duration(circleMoveDuration) // change this number to speed up or slow down the animation
        current
          .selectAll("rect.FP")
          .transition()
          .attr("x", function(d_) {
            return xScale(d_[xvar]) - symbolsize
          })
          .attr("y", function(d_) {
            return yScale(d_[yvar]) - symbolsize
          })
          .duration(circleMoveDuration) // change this number to speed up or slow down the animation
        current
          .selectAll("text.FP")
          .transition()
          .attr("x", function(d_) {
            return xScale(d_[xvar]) - symbolsize / 2
          })
          .attr("y", function(d_) {
            return yScale(d_[yvar]) + symbolsize / 2
          })
          .duration(circleMoveDuration) // change this number to speed up or slow down the animation
      })

      // these two lines cause the transition animation on the axes... they are also cause
      // chopiness in the user interface when the user slides the range sliders
      // on the right side...  uncomment to see their effect.
      svg.select(".x.axis.bottom").call(xAxisBottom)
      svg.select(".y.axis.left").call(yAxisLeft)
    }

    if (svgEnter[0][0]) {
      // dynamically generate filter sliders based on "filters" object
      $.each(dataFilters, function(i, v) {
        if (i === "agg") {
          return true
        }
        $(`<div id='${i}' class='noUi-slider'/>`).appendTo("#sliders")
        let slider = document.getElementById(i)
        $(
          `<label class='noUi-slider-label' for=${i}>${strings[i]}</label>`
        ).appendTo(slider)

        const formatttip = {
          to: function(value) {
            if (value < 1) {
              return value
            }
            if (value >= 10000) {
              return `${Math.round(value / 1000)}k`
            }
            return Math.round(value)
          },
        }

        noUiSlider.create(slider, {
          start: [v[0], v[1]], // 4 handles, starting at...
          connect: true, // Display a colored bar between the handles
          behaviour: "tap-drag", // Move handle on tap, bar is draggable
          step: v[2],
          tooltips: [formatttip, formatttip],
          range: { min: v[0], max: v[1] },
        })

        let resizeTimer
        // update filter settings when user changes slider
        slider.noUiSlider.on("update", function() {
          const filtID = this.target.id
          clearTimeout(resizeTimer)
          resizeTimer = setTimeout(function() {
            slider = document.getElementById(filtID)
            const data = slider.noUiSlider.get()
            dataFilters[filtID][0] = parseFloat(data[0])
            dataFilters[filtID][1] = parseFloat(data[1])
            plot()
          }, 25)
        })
        return true
      })

      $("#aggselect").change(function() {
        const aggchoice = $(this).val()
        dataFilters.agg = aggchoice
        plot()
      })

      $("#Xradio label").click(function() {
        currentX = $(this)
          .children("input")
          .val()
        plot()
      })
      $("#Yradio label").click(function() {
        currentY = $(this)
          .children("input")
          .val()
        plot()
      })
    }

    plot()
  }

  // easter egg
  // $("#doalittledance").click(function() {
  //   doalittledance(1600)
  // })

  chart.data = function(value) {
    if (!arguments.length) return FPdata

    value.forEach(function(d) {
      d.em_max = +d.em_max // typing these variables here for simplicity of code later on
      d.ex_max = +d.ex_max
      d.ext_coeff = +d.ext_coeff
      d.qy = +d.qy
      d.brightness = +d.brightness
    })

    FPdata = value

    saturationScale.domain([0, 35])
    return chart
  }

  // function doalittledance(_int) {
  //   const int = _int | 1000
  //   var s = ["qy", "ext_coeff", "em_max", "ex_max", "brightness"]
  //   danceInterval = setInterval(function() {
  //     var x = s[Math.floor(Math.random() * s.length)]
  //     do {
  //       var y = s[Math.floor(Math.random() * s.length)]
  //     } while (x === y)
  //     //plot(x,y);
  //   }, int)
  // }

  // function stopdancing() {
  //   clearInterval(danceInterval)
  // }

  return chart
}

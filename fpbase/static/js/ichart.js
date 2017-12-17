
///////////////////////////////////////////////////////////////
//////				Interactive Chart Logic				///////
///////////////////////////////////////////////////////////////


//  https://github.com/JacobLett/IfBreakpoint/blob/master/if-b4-breakpoint.js
// set variables
// var xs;
// var sm;
// var md;
// var lg;
// var xl;
// var breakpoint;

// // Checks if the span is set to display lock via CSS
// function checkIfBlock (target) {
//     var target = $(target).css('display') == 'block';
//     return target;
// }

// // function to check the sizes
// function checkSize (){
//   // Set some variables to use with the if checks below

// 	xs = checkIfBlock('.breakpoint-check .xs');
// 	sm = checkIfBlock('.breakpoint-check .sm');
// 	md = checkIfBlock('.breakpoint-check .md');
// 	lg = checkIfBlock('.breakpoint-check .lg');
// 	xl = checkIfBlock('.breakpoint-check .xl');


// 	// add the breakpoint to the console
// 	if( xs == true) {
// 		breakpoint = "xs - <576px";
// 		$("body").removeClass('xs sm md lg xl').addClass( "xs" );
// 	}

// 	if( sm == true) {
// 		breakpoint = "sm - ≥576px";
// 		$("body").removeClass('xs sm md lg xl').addClass( "sm" );
// 	}

// 	if( md == true) {
// 		breakpoint = "md - ≥768px";
// 		$("body").removeClass('xs sm md lg xl').addClass( "md" );
// 	}

// 	if( lg == true) {
// 		breakpoint = "lg - ≥992px";
// 		$("body").removeClass('xs sm md lg xl').addClass( "lg" );
// 	}

// 	if( xl == true) {
// 		breakpoint = "xl - ≥1200px";
// 		$("body").removeClass('xs sm md lg xl').addClass( "xl" );
// 	}

// }
// end check size


// $(document).ready(function(){
//  	// Add some invisible elements with Bootstrap CSS visibile utility classes
// 	$( "body" ).append( "<div style='display:none;' class='breakpoint-check'><span class='xs d-block d-sm-inline'></span><span class='sm d-sm-block d-md-inline'></span><span class='md d-md-block d-lg-inline'></span><span class='lg d-lg-block d-xl-inline'></span><span class='xl d-xl-block'></span></div>" );
// 	checkSize();
// });


// // Reload demo on  window resize
// $( window ).resize( function(){
// 	checkSize();
// });



//global variables to hold the current variables plotted on each axis
var currentX = "lambda_ex"
var currentY = "lambda_em"
var symbolsize = 8; //radius of circle
var bigscale = 1.5; //how much to scale up on mouseover
//global varable to set the ranges over which the data is filtered.
var filters = {
	"lambda_ex" : [350,800,1],		// array values represent [min range, max range, step (for the range slider)]
	"lambda_em" : [350,800,1],
	"E"			: [10000,170000,1000],
	"QY"		: [0,1,0.01],
	"brightness": [0,125,1]
}
//string variables for updating the axis labels
var strings = {
	"lambda_em" : "Emission Wavelength (nm)",
	"lambda_ex" : "Excitation Wavelength (nm)",
	"stokes"	: "Stokes Shift (nm)",
	"E"			: "Extinction Coefficient",
	"QY"		: "Quantum Yield",
	"brightness": "Brightness",
	"pka" 		: "pKa",
	"bleach" 	: "Bleaching Half-life (s)",
	"mature" 	: "Maturation Half-time (min)",
	"lifetime" 	: "Lifetime (ns)",
}

//shorter strings for the table
var tableStrings = {
	"Name"		: "Protein",
	"lambda_ex" : "&lambda;<sub>ex</sub> (nm)",
	"lambda_em" : "&lambda;<sub>em</sub> (nm)",
	"E"			: "EC",
	"QY"		: "QY",
	"brightness": "Brightness",
	"pka" 		: "pKa",
	"bleach" 	: "Bleaching (s)",
	"mature" 	: "Maturation (min)",
	"lifetime" 	: "Lifetime (ns)",
	"RefNum"	: "Reference"
}

// //Protein classes for tables
// var FPgroups = [
// 		{"Name" : "UV", "ex_min" : 0, "ex_max" : 380, "em_min" : 0, "em_max" : 1000, "color" : "#C080FF"},
// 		{"Name" : "Blue", "ex_min" : 380, "ex_max" : 421, "em_min" : 0, "em_max" : 470, "color" : "#8080FF"},
// 		{"Name" : "Cyan", "ex_min" : 421, "ex_max" : 473, "em_min" : 0, "em_max" : 530, "color" : "#80FFFF"},
// 		{"Name" : "Green", "ex_min" : 473, "ex_max" : 507, "em_min" : 480, "em_max" : 530, "color" : "#80FF80"},
// 		{"Name" : "Yellow", "ex_min" : 507, "ex_max" : 531, "em_min" : 500, "em_max" : 1000, "color" : "#FFFF80"},
// 		{"Name" : "Orange", "ex_min" : 531, "ex_max" : 555, "em_min" : 530, "em_max" : 569, "color" : "#FFC080"},
// 		{"Name" : "Red", "ex_min" : 555, "ex_max" : 600, "em_min" : 570, "em_max" : 620, "color" : "#FFA080"},
// 		{"Name" : "Far Red", "ex_min" : 585, "ex_max" : 631, "em_min" : 620, "em_max" : 1000, "color" : "#FF8080"},
// 		{"Name" : "Near IR", "ex_min" : 631, "ex_max" : 800, "em_min" : 661, "em_max" : 1000, "color" : "#B09090"},
// 		{"Name" : "Sapphire-type", "ex_min" : 380, "ex_max" : 420, "em_min" : 480, "em_max" : 530, "color" : "#8080FF"},
// 		{"Name" : "Long Stokes Shift", "ex_min" : 350, "ex_max" : 500, "em_min" : 560, "em_max" : 650, "color" : "#80A0FF"}
// ]

//on page load, listen to slider events and respond by updating the filter ranges (and updating the ui)
//this uses jQuery and jQuery UI which have been added to the head of the document.

$(function() {

	//dynamically generate filter sliders based on "filters" object
	$.each(filters, function(i,v){
		$("<div id='"+i+"' class='noUi-slider'/>").appendTo("#sliders");
		slider = document.getElementById(i)
		var label = $("<label class='noUi-slider-label' for="+i+">"+strings[i]+"</label>").appendTo(slider);

		noUiSlider.create(slider, {
			start: [ v[0], v[1]], // 4 handles, starting at...
			connect: true, // Display a colored bar between the handles
			behaviour: 'tap-drag', // Move handle on tap, bar is draggable
			step: v[2],
			tooltips: true,
			range: {min: v[0], max: v[1]},
			format: {
				  to: function ( value ) {
					return value;
				  },
				  from: function ( value ) {
					return value;
				  }
				}
		});

		// update filter settings when user changes slider
		slider.noUiSlider.on("update", function(){
			var filtID = this.target.id;
			slider = document.getElementById(filtID)
			data = slider.noUiSlider.get();
		 	filters[filtID][0] = parseFloat(data[0]);
		 	filters[filtID][1] = parseFloat(data[1]);
		 	plot();
		});

	});

    $( "#Xradio label").click(function() {
	  currentX = $(this).children('input').val();
	  plot();
	});
	$( "#Yradio label").click(function() {
	  currentY = $(this).children('input').val();
	  plot();
	});

	//easter egg
	$("#doalittledance").click(function(){doalittledance(1600);});
	});


// Chart dimensions.
var margin = {top: 20, right: 30, bottom: 20, left: 50},
width = 700 - margin.right,
height = 700 - margin.top - margin.bottom;

//Scales and axes
var xScale = d3.scale.linear()
			.range ([0, width]);

var yScale = d3.scale.linear()
			.range ([height, 0]);

//This scale will set the saturation (gray to saturated color).  We will use it for mapping brightness.
var saturationScale = d3.scale.linear()
			.range([0, 1])
			.domain([0, 100]);

//This scale will set the hue.  We will use it for mapping emission wavelength.
var hueScale = d3.scale.linear()
			.range([300, 300, 240, 0, 0])
			.domain([200, 405, 440, 650, 850]);

//X and Y axes
var xAxis_bottom = d3.svg.axis().scale(xScale).tickSize(5).tickSubdivide(true);
var yAxis_left = d3.svg.axis().scale(yScale).tickSize(5).orient("left").tickSubdivide(true);

//top and right axes are identical but without tick labels
var xAxis_top = d3.svg.axis().scale(xScale).tickSize(5).orient("top").tickSubdivide(true).tickFormat(function (d) { return ''; });;;
var yAxis_right = d3.svg.axis().scale(yScale).tickSize(5).orient("right").tickSubdivide(true).tickFormat(function (d) { return ''; });;

// Create the SVG container and set the origin.
var svg = d3.select("#graph").append("svg")
	//.attr("width", width + margin.left + margin.right)
	//.attr("height", height + margin.top + margin.bottom)
	.attr("viewBox", "0 0 760 760")
	.attr("id", "mainchart")
	.attr("preserveAspectRatio", "xMinYMin meet")
	.append("g")
	.attr("transform", "translate(" + margin.left + "," + margin.top + ")")
	.classed("svg-content-responsive", true);

//Add the axes
svg.append("g")
	.attr("class", "x axis bottom")
	.attr("transform", "translate(0," + height + ")")
	.call(xAxis_bottom);
svg.append("g")
	.attr("class", "y axis left")
	.call(yAxis_left);
svg.append("g")
	.attr("class", "x axis top")
	.call(xAxis_top);
svg.append("svg:g")
	.attr("class", "y axis right")
	.attr("transform", "translate(" + width + ",0)")
	.call(yAxis_right);

// Add an x-axis label.
svg.append("text")
	.attr("class", "x label")
	.attr("text-anchor", "middle")
	.attr("x", width/2 )
	.attr("y", height-10)
	.text("Excitation wavelength (nm)");

// Add a y-axis label.
svg.append("text")
	.attr("class", "y label")
	.attr("text-anchor", "middle")
	.attr("x", -height/2)
	.attr("y", margin.left-30)
	.attr("transform", "rotate(-90)")
	.text("Emission wavelength (nm)");

//Add a clipping path so that data points don't go outside of frame
svg.append("clipPath")                  //Make a new clipPath
	.attr("id", "chart-area")           //Assign an ID
		.append("rect")
		.attr("width", width)
		.attr("height", height);

//enable zooming
var zoom = d3.behavior.zoom()
	.x(xScale)
	.y(yScale)
	.scaleExtent([1, 10])
	.on("zoom", draw_graph);

function plotcircle(sel){
	circle = sel.append("circle")
		.attr("class", "FP")
		.attr("r", symbolsize)
		.attr("stroke", "#000")
		.attr("opacity", 0.7)
		.style("fill", function (d) { return d3.hsl(hueScale (d.lambda_em), saturationScale (d.brightness), 0.5)});
		addactions(circle);
	}

function plotsquare(sel){
	square = sel.append("rect")
		.attr("class", "FP")
		.attr("width", symbolsize*2)
		.attr("height", symbolsize*2)
		.attr("stroke", "#000")
		.attr("opacity", 0.7)
		.style("fill", function (d) { return d3.hsl(hueScale (d.lambda_em), saturationScale (d.brightness), 0.5)});
		addactions(square);
	}

	function plottext(sel){
	text = sel.append("text")
		.attr("class", "FP")
		.text(function (d) {
			if (d["agg"] == "d") { return "2"}
			else if (d["agg"] == "td") { return "t"}
			else if (d["agg"] == "t") { return "4"}
		;} )
	}

function addactions(sel){
		sel.on('click', function(e){
			window.location = e.url;
		})
		sel.on("mouseover", function(d) {
			//Get this bar's x/y values, then augment for the tooltip
			if (d3.select(this).attr("cx")){ //if circle
				d3.select(this).transition().duration(100).attr("r",symbolsize*bigscale);
				var xPosition = parseFloat(d3.select(this).attr("cx"))
				var yPosition = parseFloat(d3.select(this).attr("cy"))
			} else if (d3.select(this).attr("x")){ //if rectangle
				d3.select(this).transition().duration(100)
					.attr("x", function (d) { return xScale (d[currentX]) - symbolsize*bigscale; })
					.attr("y", function (d) { return yScale (d[currentY]) - symbolsize*bigscale; })
					.attr("width", symbolsize*2*bigscale)
					.attr("height", symbolsize*2*bigscale);
				var xPosition = parseFloat(d3.select(this).attr("x") )
				var yPosition = parseFloat(d3.select(this).attr("y") )
			}
			if (xPosition<width*2/3){
				xPosition += 70;
			} else {
				xPosition -= 140;
			}
			if (yPosition>520){
				yPosition = 520;
			}
			//Update the tooltip position and value
			d3.select("#tooltip")
				.style("left", xPosition + "px")
				.style("top", yPosition + "px")
				.select("#exvalue")
				.text(d.lambda_ex)
			d3.select("#tooltip")
				.select("#emvalue")
				.text(d.lambda_em);
			d3.select("#tooltip")
				.select("#ecvalue")
				.text(d.E);
			d3.select("#tooltip")
				.select("#qyvalue")
				.text(d.QY);
			d3.select("#tooltip")
				.select("h3")
				.html(d.name);
			d3.select("#tooltip")
				.select("#brightnessvalue")
				.text(d.brightness);

		//Show the tooltip
		d3.select("#tooltip").classed("hidden", false);
		})

		.on("mouseout", function() {
			if (d3.select(this).attr("cx")){ //if circle
				d3.select(this).transition().duration(200).attr("r",symbolsize)
			} else if (d3.select(this).attr("x")){ //if circle
				d3.select(this).transition().duration(200)
					.attr("x", function (d) { return xScale (d[currentX]) - symbolsize; })
					.attr("y", function (d) { return yScale (d[currentY]) - symbolsize; })
					.attr("width", symbolsize*2)
					.attr("height", symbolsize*2);
			}
			//Hide the tooltip
			d3.select("#tooltip").classed("hidden", true);
		})
		}

svg.append("rect")
	.attr("class", "pane")
	.attr("width", width)
	.attr("height", height)
	.call(zoom);

var FPdata = []; //Where the fluorescent protein data table will end up.

function draw_graph(){
	//redraw axes with new domains
	svg.select(".x.axis.bottom").call(xAxis_bottom);
	svg.select(".y.axis.left").call(yAxis_left);
	svg.select(".x.axis.top").call(xAxis_top);
	svg.select(".y.axis.right").call(yAxis_right);

	svg.selectAll("circle.FP")
		.attr("cx", function (d) { return xScale (d[currentX]); })
		.attr("cy", function (d) { return yScale (d[currentY]); })

	svg.selectAll("rect.FP")
	    .attr("x", function (d) { return xScale (d[currentX]) - symbolsize; })
	    .attr("y", function (d) { return yScale (d[currentY]) - symbolsize; })

	svg.selectAll("text.FP")
	    .attr("x", function (d) { return xScale (d[currentX]) - symbolsize/2; })
	    .attr("y", function (d) { return yScale (d[currentY]) + symbolsize/2; })
}

//i added this more flexible plotting function to be able to plot different variables on each axis.  It takes three optional parameters: the data array, and two axes variables.
function plot(xvar,yvar,data){
	//set default values... if plot() is called without arguments, these default values will be used.
	xvar = xvar || currentX;
	yvar = yvar || currentY;
	data = data || FPdata;

	// helper function to iterate through all of the data filters (without having to type them all out)
	function filtercheck(data){
		for (f in filters){
			v = filters[f];
			if( data[f] < v[0] || data[f] > v[1] ) {return false;}
		}
		return true;
	}

	//filter the data according to the user settings for EC, QY, and brightness range
	data = data.filter(function(d) {  return filtercheck(d) ? d : null; });

	//filter out data with empty values
	data = data.filter(function(d) {return d[xvar] > 0 && d[yvar] > 0;});

	//update scale domains based on data
	xScale.domain([
		d3.min (data, function(d) { return .99 * d[xvar]; }),
		d3.max (data, function(d) { return 1.01 * d[xvar]; })
	])
	.nice();
	zoom.x(xScale);

	yScale.domain([
		d3.min (data, function(d) { return .99 * d[yvar]; }),
		d3.max (data, function(d) { return 1.01 * d[yvar]; })
	])
	.nice();
	zoom.y(yScale);

	//relabel X and Y axes
	svg.select(".x.label").text(strings[xvar])
	svg.select(".y.label").text(strings[yvar])

	// Join new data with old elements, if any.
	var datagroup = svg.selectAll("g.FP").data(data, function (d){ return d.name;});
	entergroup = datagroup.enter().append("g")
		.attr("class", "FP")
		.attr("clip-path", "url(#chart-area)")
		.call(zoom);		//so we can zoom while moused over elements

	entergroup.each(function(d, i) {
		//determine type of protein and whether to plot a circle or a square
		if (d["type"] =="e"){
			//plot squeares (for proteins with cofactor)
			plotsquare(d3.select(this));
		} else {
			// plot new squares
			plotcircle(d3.select(this));
		}
		//add text to markers
		plottext(d3.select(this));
	})

	// Remove old elements as needed.
	datagroup.exit().remove();

	// move circles to their new positions (based on axes) with transition animation
	datagroup.each(function(d, i) {
		current = d3.select(this)
		current.selectAll("circle.FP")
			.transition()
			.attr("cx", function (d) { return xScale (d[xvar]); })
			.attr("cy", function (d) { return yScale (d[yvar]); })
			.duration(800); //change this number to speed up or slow down the animation
		current.selectAll("rect.FP")
			.transition()
			.attr("x", function (d) { return xScale (d[xvar]) - symbolsize; })
			.attr("y", function (d) { return yScale (d[yvar]) - symbolsize; })
			.duration(800); //change this number to speed up or slow down the animation
		current.selectAll("text.FP")
			.transition()
			.attr("x", function (d) { return xScale (d[xvar]) - symbolsize/2; })
			.attr("y", function (d) { return yScale (d[yvar]) + symbolsize/2; })
			.duration(800); //change this number to speed up or slow down the animation
	})

	// these two lines cause the transition animation on the axes... they are also cause chopiness in the user interface when the user slides the range sliders on the right side...  uncomment to see their effect.
	svg.select(".x.axis.bottom").call(xAxis_bottom);
	svg.select(".y.axis.left").call(yAxis_left);
}


function doalittledance(int) {
	var s = ["QY","E","lambda_em","lambda_ex","brightness"];
	setInterval(function() {
	  var x = s[Math.floor(Math.random() * s.length)];
	  do{
	    var y = s[Math.floor(Math.random() * s.length)];
	  }	while (x == y);
	  plot(x,y);
	}, int);

}


// //this bit is just a jQuery plugin to make the radio checkboxes on the right side vertical
// (function( $ ){
// //plugin buttonset vertical
// $.fn.buttonsetv = function() {
//   $(':radio, :checkbox', this).wrap('<div style="margin: 1px"/>');
//   $(this).buttonset();
//   $('label:first', this).removeClass('ui-corner-left').addClass('ui-corner-top');
//   $('label:last', this).removeClass('ui-corner-right').addClass('ui-corner-bottom');
//   mw = 0; // max witdh
//   $('label', this).each(function(index){
//      w = $(this).width();
//      if (w > mw) mw = w;
//   })
//   $('label', this).each(function(index){
//     $(this).width(mw);
//   })
// };
// })( jQuery );


///////////////////////////////////////////////////////////////
//////				END Interactive Chart Logic			///////
///////////////////////////////////////////////////////////////
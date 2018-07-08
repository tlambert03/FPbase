function createDropdown(data){

	var dropdown = d3.select("body").append("select")
			.on("change", function(){

				var selection = dropdown.property('selectedIndex');

				myData =[]
				myData.push({
					key 	: data[selection].name + ' ex',
					values 	: formatForNVD3(data[selection].ex_spectrum),
					peak	: data[selection].ex_max,
					type 	: 'ex',
					color 	: realColor(data[selection].ex_max),
					area 	: true
				});
				myData.push({
					key 	: data[selection].name + ' em',
					values 	: formatForNVD3(data[selection].em_spectrum),
					peak	: data[selection].em_max,
					type 	: 'em',
					color 	: realColor(data[selection].em_max),
					area 	: true
				});


				svg.datum(myData)
				chart.update();
			});

	list = dropdown.selectAll('option').data(data);

	list.enter()
		.append("option")
		.attr("value", function(d){return d.name})
		.text(function(d){return d.name});

};

function createRadioBoxes(data){

	var radios = d3.select("#protein_form").append("form");

	list = radios.selectAll('label').data(data);

	var listitem = list.enter()
		.append('label')
		    .attr('for', function(d){return d.name});

	listitem.append('input')
			.attr('type', 'checkbox')
			.attr('name', 'protein')
			.attr('value', function(d){return d.name})
			.attr('id', function(d){return d.name})
			.on("click", function(d){
				myData =[]
				data = d3.selectAll('input:checked').data();

				var arrayLength = data.length;
				for (var i = 0; i < arrayLength; i++) {
				    myData.push({
						key 	: data[i].name + ' ex',
						values 	: formatForNVD3(data[i].ex_spectrum),
						peak	: data[i].ex_max,
						type 	: 'ex',
						color 	: realColor(data[i].ex_max),
						area 	: true
					});
					myData.push({
						key 	: data[i].name + ' em',
						values 	: formatForNVD3(data[i].em_spectrum),
						peak	: data[i].em_max,
						type 	: 'em',
						color 	: realColor(data[i].em_max),
						area 	: true
					});
				}

				svg.datum(myData)
				chart.update();

			});

	listitem.append('text')
			.text(function(d){return d.name})
			.append('span').html('<br>');
}

function formatForNVD3(input){

	output = [];
	var arrayLength = input.length;
	var first = input[0][0];
	if (first>300){
		for(n=0; n < (first-300) ; n++){
			output.push({x: (300+n) , y: 0})
		}
	}
	for (var i = 0; i < arrayLength; i++) {
	    output.push({x: input[i][0],y: input[i][1]})
	    //Do something
	}
	var leftover = 800-input[arrayLength-1][0]

	if (leftover>0){
		for(n=0; n < leftover+1 ;n++ ){
			output.push({x: (800-leftover+n), y: 0})
		}
	}

	return output;
}

function sortBy(prop,withm) {
    return function(a, b) {
    		if (prop == 'name') {
    			var aComp = a[prop]
    			    bComp = b[prop]
    			// by default, ignore the "m" at the beginning of proteins when sorting unless the sort is called using "withm=true"
    			if (!withm) {
    				if (a[prop][0] == 'm') {aComp = aComp.slice(1)}
    				if (b[prop][0] == 'm') {bComp = bComp.slice(1)}
    			}
    			if (aComp.toLowerCase() < bComp.toLowerCase()) return -1;
    			if (aComp.toLowerCase() > bComp.toLowerCase()) return 1;
    		}

		return a[prop] - b[prop];
    }
}

function SpectraFileTypes(d){
  return {
	name: d.fluor_name,
	ex_spectra: JSON.parse(d.ex_spectrum),
	em_spectra: JSON.parse(d.em_spectrum),
	ex_max: +d.ex_max,
	em_max: +d.em_max,
  }
}

var hueScale = d3.scale.linear()
				.domain([200, 380, 500, 640, 850])
				.range([300, 300, 180,0, 0]);


function realColor(d) {
	return d3.hsl(hueScale(d), 1, 0.5);
}

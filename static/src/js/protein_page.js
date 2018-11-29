import LineageChart from './lineage.js';

window.initSimpleSpectra = function(selection, myData, protein) {
  // disable 2p by default
  for (var n = 0; n < myData.length; n++) {
    if (myData[n].key.indexOf('2p') >= 0 && myData.length > 1) {
      myData[n].disabled = true;
    }
  }

  var interactive = true;
  if (typeof window.mobilecheck === 'function') {
    interactive = !window.mobilecheck();
  }
  var chart = nv.models.lineChart();
  chart.yDomain([0, 1.04]);
  var svg = d3.select(selection);

  if (protein){
    $(selection).prepend(
      $('<desc>', { id: 'svgDesc' }).text(
        `Fluorescent protein ${protein} excitation and emission spectra`
      )
    );
    $(selection).prepend(
      $('<title>', { id: 'svgTitle' }).text(`${protein} Spectrum`)
    );
  }

  chart.options({
    margin: { left: 20, bottom: 30 },
    noData: 'Unable to load data...',
    transitionDuration: 200,
    showLegend:
      myData.length > 2 ||
      myData.some(function(a) {
        return a.type == '2p';
      }),
    showXAxis: true,
    showYAxis: false,
    useInteractiveGuideline: interactive,
    forceY: [0, 1.04],
    forceX: [350, 750],
  });

  /*These lines are all chart setup.  Pick and choose which chart features you want to utilize. */
  nv.addGraph(
    function() {
      chart.xAxis //Chart x-axis settings
        //.axisLabel('Wavelength (nm)')
        .tickFormat(d3.format(',r'));

      // chart.yAxis     //Chart y-axis settings
      //   .axisLabel('Normalized Ex/Em')
      //   .tickFormat(d3.format('1%'));

      svg
        .datum(myData) //Populate the <svg> element with chart data...
        .call(chart); //Finally, render the chart!

      //Update the chart when window resizes.

      nv.utils.windowResize(function() {
        chart.update();
      });
      return chart;
    },
    function() {
      setTimeout(function() {
        chart.update();
      }, 50);
    }
  );
};

window.initSnapGene = function(protein, selection) {
  var snaptemplate =
    'https://www.snapgene.com/resources/plasmid_files/fluorescent_protein_genes_and_plasmids/*/';
  $.get({
    url: snaptemplate.replace('*/', 'plasmids.xml'),
    success: function(data) {
      /* handle data here */
      var hits = [];
      $(data)
        .find(`[name*="${protein}"]`)
        .each(function() {
          var name = $(this).attr('name');
          while (name.charAt(0) === 'p') {
            name = name.substr(1);
          }
          $.each(
            ['-N2', '-N1', '-N', '-C2', '-C1', '-C', '-1', '-B', '-S1', '-PRL'],
            function(index, val) {
              name = name.replace(val, '');
            }
          );
          if (name == protein || name == `${protein}1`) {
            hits.push([$(this).attr('name'), $(this).attr('url')]);
          }
        });

      if (hits.length) {
        var t = 'SnapGene links: ';
        if (hits.length == 1) {
          t = 'SnapGene link: ';
        }
        var $li = $('<li>')
          .text(t)
          .appendTo($(selection));
        $.each(hits, function(index, val) {
          $li.append(
            $('<a>', { href: snaptemplate.replace('*', val[1]) }).text(val[0])
          );
          if (index != hits.length - 1) {
            $li.append(', ');
          }
        });
      }
    }
  });
};

$(window).on('load', function() {
  $('#excerptModal').on('show.bs.modal', function () {
      $('#excerptModal').css('z-index', 1200);
  });

  $('#excerptModal').on('hidden.bs.modal', function () {
      $('#excerptModal').css('z-index', 1039);
  });

  var resizeTimer;
  $('.lineage').each(function(i, el) {
    var slug = $(el).data('slug') || '';
    $.getJSON('/ajax/lineage/' + slug, function(data) {
      var linchart = LineageChart({slug: slug}).data(data);
      var lineage = d3.select(el);
      lineage.call(linchart);
      $(window).on('resize', function(e) {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(function() {
          lineage.call(linchart);
        }, 250);
      });
      linchart.duration(200);
      if (slug && $('.lineage-wrapper')){
        var slugpos = d3.select(`#${slug}_circle`).datum().y;
        $('.lineage-wrapper').scrollLeft(slugpos - ((window.innerWidth - 30) / 3));
      }

    });
  });
});

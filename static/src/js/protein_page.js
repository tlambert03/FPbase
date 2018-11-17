import LineageChart from './lineage.js'
import initPDB from './my-litemol.js'

window.initPDB = initPDB;

$(function(){

  var resizeTimer;
  $('.lineage').each(function(i, el){
    var slug = $(el).data('slug') || '';
    $.getJSON("/ajax/lineage/" + slug, function(data){
      var linchart = LineageChart().data(data);
      var lineage = d3.select(el)
      lineage.call(linchart, slug);
      $(window).on('resize', function(e) {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(function() {
          lineage.call(linchart);
        }, 250);
      });
    })
  });

});

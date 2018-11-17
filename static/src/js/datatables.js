import 'datatables.net-bs4';


///////////DATA TABLE

$(function() {
  $('#proteinTable').DataTable({
    "pageLength": 50,
    "autoWidth": false,
    "lengthMenu": [ [10, 25, 50, 100, -1], [10, 25, 50, 100,"All"] ],
    'fixedHeader': {
        'header': false,
        'footer': true
    },
    "order": [[ 6, 'desc' ]],
  });

  $('.table-filter').change(function(){
    var searchval = this.value;
    if (searchval != ''){
      searchval = '^' + this.value +'$';
    }
    searchcol = $(this).attr('data-col')
    $('#proteinTable').DataTable().column('.' + searchcol).search(searchval, true, false).draw();
  });

  // $('.exmax').each(function(){
  //  $(this).closest('tr').css('background-color', get_color_group($(this).html()))
  // })

});

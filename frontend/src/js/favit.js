$(document).ready(function() {
  $('#add_remove_favorite').click(function(e) {
      var $obj = $(this);
      var target_id = $obj.data('target').split('_')[1];
      $obj.prop('disabled', true);
      $.ajax({
      url: $obj.attr('data-action-url'),
      type: 'POST',
      data: {target_model: $obj.data('model'),
             target_object_id: $obj.data('target').split('_')[1],
             csrfmiddlewaretoken: window.CSRF_TOKEN},
      success: function(response) {
          if (response.status == 'added') {
            $obj.find('svg').attr('data-prefix', 'fas');
            $obj.find('span').text('Remove from favorites');
          }
          else {
            $obj.find('svg').attr('data-prefix', 'far');
            $obj.find('span').text('Add to favorites');
          }
          $obj.parent('.favit').children('.fav-count').text(response.fav_count);
          $obj.prop('disabled', false);
      }
      });
      e.preventDefault(); // avoid to execute the actual submit of the form.
  });


  $('.btn.unfave').click(function() {
    var $obj = $(this);
    $obj.prop('disabled', true);
    $.ajax({
      url: $obj.attr('data-action-url'),
      type: 'POST',
      data: {
        target_model: $obj.data('model'),
        target_object_id: $obj.data('id'),
        csrfmiddlewaretoken:  window.CSRF_TOKEN,
      },
      success: function(response) {
        if (response.status == 'deleted') {
          $obj.parent().remove();
        }
      },
      complete: function(response) {
        $obj.prop('disabled', false);
      }
    });
    e.preventDefault(); // avoid to execute the actual submit of the form.
  });
});

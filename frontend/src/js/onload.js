$(function() {
    $("#protein-image .no-glow").fadeOut(2200)
    $("#protein-image .glow").fadeIn(2200)

    $('[data-toggle="tooltip"]').tooltip({
      trigger : 'hover',
      delay: { "show": 200 },
    });

    $('[data-toggle="popover"]').popover({ html: true });
    $('.popover-dismiss').popover({ trigger: 'focus' });

});

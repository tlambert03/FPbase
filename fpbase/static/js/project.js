/* Project specific Javascript goes here. */

/*
Formatting hack to get around crispy-forms unfortunate hardcoding
in helpers.FormHelper:

    if template_pack == 'bootstrap4':
        grid_colum_matcher = re.compile('\w*col-(xs|sm|md|lg|xl)-\d+\w*')
        using_grid_layout = (grid_colum_matcher.match(self.label_class) or
                             grid_colum_matcher.match(self.field_class))
        if using_grid_layout:
            items['using_grid_layout'] = True

Issues with the above approach:

1. Fragile: Assumes Bootstrap 4's API doesn't change (it does)
2. Unforgiving: Doesn't allow for any variation in template design
3. Really Unforgiving: No way to override this behavior
4. Undocumented: No mention in the documentation, or it's too hard for me to find
*/
$('.form-group').removeClass('row');

//Navbar Scroll Event.
// comment this out to remove "navbar hiding" when the user scrolls down
// var lastScrollTop = 0;
// var navbar        = $('.navbar');
// $(window).scroll(function(event){
//    var st = $(this).scrollTop();
//    if (st > lastScrollTop){
//        navbar.addClass('navbar-scroll-custom');
//    } else {
//       navbar.removeClass('navbar-scroll-custom');
//    }
//    lastScrollTop = st;
// });


///////////////////////
// FORM BEHAVIOR
///////////////////////

// Hide unnecessary fiels in forms when the dark state checkbox is toggled
$('.dark_state_button input').click(function(){
	neighbors = $(this).closest('.stateform_block').children('.hide_if_dark')
	if ($(this)[0].checked){
		neighbors.hide()
        $(this).find("input[name*='max']").empty()
	}
	else{
		neighbors.show()
        $(this).closest('.stateform_block')
	}
});

$(function() {
    $('.dark_state_button input:checked').each(function(){
        neighbors = $(this).closest('.stateform_block').children('.hide_if_dark')
        neighbors.hide()
    });
});

function reset_ipgid(hintstring){
    if (hintstring){
        $("#hint_id_ipg_id").html(hintstring)
    }
    $("#id_seq").prop('disabled', false);
    $("#hint_id_seq").html('Amino acid sequence (IPG ID is preferred)')
}

$("#id_ipg_id").change(function(){
    ipg_id = $(this).val()
    protein_uri = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=json&rettype=fasta&id="
    ipg_uri = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=ipg&retmode=json&id="
    fpbase_params = "&tool=fpbase&email=talley.lambert+fpbase@gmail.com"
    $.ajax({
        url: ipg_uri + ipg_id + fpbase_params,
        context: document.body,
        success: function(data){
            if (!('result' in data)){
                reset_ipgid('NCBI <a href="https://www.ncbi.nlm.nih.gov/ipg/docs/about/">Identical Protein Group ID</a>')
            } else if (ipg_id in data.result){
                var accession = data.result[ipg_id].accession
                title = data.result[ipg_id].title
                $("#hint_id_ipg_id").html("IPG name: " + title)
                $.ajax({
                    url: protein_uri + accession + fpbase_params,
                    context: document.body,
                    success: function(data2){
                        var lines = data2.split('\n');
                        var seq = ''
                        for(var i = 0;i < lines.length;i++){
                            if (lines[i].length != 0 && lines[i][0] != '>') {
                                seq += lines[i]
                            }
                        }
                        $("#id_seq").val(seq)
                        $("#id_seq").prop('disabled', true);
                        $("#hint_id_seq").html('Sequence input disabled when IPG ID provided')
                    },
                    error: function(data){
                        reset_ipgid('Unrecognized IPG ID')
                    }
                })
            }
        },
        error: function(data){
            reset_ipgid('Unrecognized IPG ID')
        }
    });

})

// auto populate PMID after DOI input
$('#id_reference_doi').change(function(){
	doi = $(this).val()

	searchurl = "https://api.crossref.org/v1/works/"
	searchurl += doi

	// searchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="
	// searchurl += doi + "[doi]&retmode=json"

	$.ajax({
	    url: searchurl,
	    context: document.body,
	    success: function(data){
	    	if (data.status == 'ok'){
	    		year = data.message.issued["date-parts"]["0"]["0"];
	    		author = data.message.author["0"]["family"];
	    		title = data.message.title[0].slice(0, 45);
	    		citation = author + " (" + year + ") " + title + "...";
	    		$("#hint_id_reference_doi").html(citation)
	    	}
	    	else{
	    		$("#hint_id_reference_doi").html("DOI not found at Crossref")
	    	}
	    },
	    error: function(data){
	    	$("#hint_id_reference_doi").html("DOI not found at Crossref")
	    }
	});
})

// // auto populate DOI after PMID input
// $('#id_reference_pmid').change(function(){
// 	pmid = $(this).val();
// 	console.log(pmid);
// 	searchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=";
// 	searchurl += pmid + "&retmode=json";
// 	$.ajax({
// 	    url: searchurl,
// 	    context: document.body,
// 	    success: function(data){
// 	    	console.log(data)
// 	    	if (pmid in data.result){
// 	    		articleids = data.result[pmid].articleids;
// 	    		doi = articleids.find(o => o.idtype === 'doi');
// 	    		if(doi){
// 	    			$('#id_reference_doi').val(doi.value);
// 	    		} else {
// 	    			$('#id_reference_doi').val('');
// 	    			$("#hint_id_reference_doi").html("DOI not found at Pubmed")
// 	    		}
// 	    		year = data.result[pmid].pubdate.split(" ")[0];
// 	    		author = data.result[pmid].sortfirstauthor.split(" ")[0];
// 	    		title = data.result[pmid].title.slice(0, 38);
// 	    		citation = author + " (" + year + ") " + title + "...";
// 	    		$("#hint_id_reference_doi").html(citation)
// 	    		$("#hint_id_reference_pmid").html("Pubmed ID located!")
// 	    	} else {
// 	    		$("#hint_id_reference_pmid").html("Unable to locate Pubmed ID")
// 	    		$('#id_reference_doi').val('');
// 	    		$("#hint_id_reference_doi").html("DOI not found at Pubmed")
// 	    	}
// 	    }
// 	});
// })



/////////////////// PROTEIN DETAIL PAGE //////////////////////


function chunkString(str, len) {
    var _size = Math.ceil(str.length/len),
        _ret  = new Array(_size),
        _offset
    ;

    for(var _i=0; _i<_size; _i++) {
      _offset = _i * len;
      _ret[_i] = str.substring(_offset, _offset + len);
    }

    return _ret;
}

function formatAAseq(elem, breakpoint) {
    breakpoint = breakpoint || 10;
    //clear any existing counts
    elem.find('.sequence_count').empty()
    // extract the string and chop it up into segments by breakpoint
    words = chunkString(elem.text().replace(/ /g,''), breakpoint)
    // clear the div
    elem.empty()
    //console.log(words)
    // create new divs for the formatted content
    seqcount = $("<div class='sequence_count'></div>").appendTo(elem);
    seqdiv = $("<div class='formatted_aminosquence'></div>").appendTo(elem);
    seqdiv.text(words[0]);
    seqcount.append(1 + '<br>')
    var height = seqdiv.height();
    for (var i = 1; i < words.length; i++) {
        seqdiv.text(seqdiv.text() + " " + words[i]);
        if (seqdiv.height() > height) {
            // line break occured at iteration i
            //console.log(words[i])
            seqcount.append((i*breakpoint) + 1 + '<br>')
            height = seqdiv.height();
        }
    }
    elem.show()
}

$(function(){
	$('.aminosequence').each(function(){
		formatAAseq($(this));
	})
})

$(window).resize(function() {
	$('.aminosequence').each(function(){
		formatAAseq($(this));
	})
});


  $("#refModalForm").submit(function(e) {
      var form = $(this).closest("form");
      $.ajax({
             type: "POST",
             url: form.attr("data-action-url"),
             data: form.serialize(),
             dataType: 'json',
             success: function(data) {location.reload();}
           });
      e.preventDefault(); // avoid to execute the actual submit of the form.
  });

  function register_transition_form(){
    $('.trans_formset_div').formset({
      addText: 'Add Transition',
      addCssClass: 'btn btn-info mb-4',
      deleteCssClass: 'transDelete',
      deleteText: '<i class="fas fa-minus-circle"></i>',
      prefix: 'transitions',
      processHidden: true,  // I added this to
    });
  }

  // This function is for showing the modal
  $(function () {
       $("#show_transition_modal").click(function () {
           $.ajax({
               type: 'GET',
               url:  $(this).attr('data-action-url'),
               data: {},
               cache: false,
               success: function (data, status) {
                   $("#transitionForm").html(data)
                   register_transition_form();
                   $('#transitionModal').modal()

               }
           });
       });
  });

  $("#transitionForm").submit(function(e) {
    var form = $(this).closest("form");
    $.ajax({
            type: "POST",
      url: form.attr("data-action-url"),
            data: form.serialize(),
            cache: false,
            success: function(data, status) {location.reload();},
            error: function(data, status, error) {
              $("#transitionForm").html(data.responseText)
              register_transition_form();
            }
    });
      e.preventDefault(); // avoid to execute the actual submit of the form.
  });


  $("#adminApprove, #adminRevert").submit(function(e) {
    var form = $(this).closest("form");
    $.ajax({
        type: "POST",
        url: form.attr("data-action-url"),
        data: form.serialize(),
        dataType: 'json',
        success: function(data) {window.location = '{{ protein.get_absolute_url }}'},
    });
      e.preventDefault(); // avoid to execute the actual submit of the form.
  });

  /////////////////// END PROTEIN DETAIL PAGE //////////////////////


window.mobilecheck = function() {
  var check = false;
  (function(a){if(/(android|bb\d+|meego).+mobile|avantgo|bada\/|blackberry|blazer|compal|elaine|fennec|hiptop|iemobile|ip(hone|od)|iris|kindle|lge |maemo|midp|mmp|mobile.+firefox|netfront|opera m(ob|in)i|palm( os)?|phone|p(ixi|re)\/|plucker|pocket|psp|series(4|6)0|symbian|treo|up\.(browser|link)|vodafone|wap|windows ce|xda|xiino/i.test(a)||/1207|6310|6590|3gso|4thp|50[1-6]i|770s|802s|a wa|abac|ac(er|oo|s\-)|ai(ko|rn)|al(av|ca|co)|amoi|an(ex|ny|yw)|aptu|ar(ch|go)|as(te|us)|attw|au(di|\-m|r |s )|avan|be(ck|ll|nq)|bi(lb|rd)|bl(ac|az)|br(e|v)w|bumb|bw\-(n|u)|c55\/|capi|ccwa|cdm\-|cell|chtm|cldc|cmd\-|co(mp|nd)|craw|da(it|ll|ng)|dbte|dc\-s|devi|dica|dmob|do(c|p)o|ds(12|\-d)|el(49|ai)|em(l2|ul)|er(ic|k0)|esl8|ez([4-7]0|os|wa|ze)|fetc|fly(\-|_)|g1 u|g560|gene|gf\-5|g\-mo|go(\.w|od)|gr(ad|un)|haie|hcit|hd\-(m|p|t)|hei\-|hi(pt|ta)|hp( i|ip)|hs\-c|ht(c(\-| |_|a|g|p|s|t)|tp)|hu(aw|tc)|i\-(20|go|ma)|i230|iac( |\-|\/)|ibro|idea|ig01|ikom|im1k|inno|ipaq|iris|ja(t|v)a|jbro|jemu|jigs|kddi|keji|kgt( |\/)|klon|kpt |kwc\-|kyo(c|k)|le(no|xi)|lg( g|\/(k|l|u)|50|54|\-[a-w])|libw|lynx|m1\-w|m3ga|m50\/|ma(te|ui|xo)|mc(01|21|ca)|m\-cr|me(rc|ri)|mi(o8|oa|ts)|mmef|mo(01|02|bi|de|do|t(\-| |o|v)|zz)|mt(50|p1|v )|mwbp|mywa|n10[0-2]|n20[2-3]|n30(0|2)|n50(0|2|5)|n7(0(0|1)|10)|ne((c|m)\-|on|tf|wf|wg|wt)|nok(6|i)|nzph|o2im|op(ti|wv)|oran|owg1|p800|pan(a|d|t)|pdxg|pg(13|\-([1-8]|c))|phil|pire|pl(ay|uc)|pn\-2|po(ck|rt|se)|prox|psio|pt\-g|qa\-a|qc(07|12|21|32|60|\-[2-7]|i\-)|qtek|r380|r600|raks|rim9|ro(ve|zo)|s55\/|sa(ge|ma|mm|ms|ny|va)|sc(01|h\-|oo|p\-)|sdk\/|se(c(\-|0|1)|47|mc|nd|ri)|sgh\-|shar|sie(\-|m)|sk\-0|sl(45|id)|sm(al|ar|b3|it|t5)|so(ft|ny)|sp(01|h\-|v\-|v )|sy(01|mb)|t2(18|50)|t6(00|10|18)|ta(gt|lk)|tcl\-|tdg\-|tel(i|m)|tim\-|t\-mo|to(pl|sh)|ts(70|m\-|m3|m5)|tx\-9|up(\.b|g1|si)|utst|v400|v750|veri|vi(rg|te)|vk(40|5[0-3]|\-v)|vm40|voda|vulc|vx(52|53|60|61|70|80|81|83|85|98)|w3c(\-| )|webc|whit|wi(g |nc|nw)|wmlb|wonu|x700|yas\-|your|zeto|zte\-/i.test(a.substr(0,4))) check = true;})(navigator.userAgent||navigator.vendor||window.opera);
  return check;
};

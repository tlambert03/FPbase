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
	}
	else{
		neighbors.show()
	}
})

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

function format_protein(elem, breakpoint) {
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
            console.log("height: " + height, "length: ", seqdiv.text().length)
        }
    }
    elem.show()
}

$(function(){
	$('.aminosequence').each(function(){
		format_protein($(this));
	})
})

$(window).resize(function() {
	$('.aminosequence').each(function(){
		format_protein($(this));
	})
});
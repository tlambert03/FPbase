(function (global) {
  // add array index of for old browsers (IE<9)
  if (!Array.prototype.indexOf) {
    Array.prototype.indexOf = function(obj, start) {
      var i, j;
      i = start || 0;
      j = this.length;
      while (i < j) {
        if (this[i] === obj) {
          return i;
        }
        i++;
      }
      return -1;
    };
  }

  // make a global object to store stuff in
  if(!global.FPbase) { global.FPbase = {}; };
  var FPbase = global.FPbase;

  // To keep track of which embeds we have already processed
  if(!FPbase.processedScripts) { FPbase.processedScripts = []; };
  var processedScripts = FPbase.processedScripts;

  if(!FPbase.styleTags) { FPbase.styleTags = []; };
  var styleTags = FPbase.styleTags;

  var scriptTags = document.getElementsByTagName('script');
  var thisRequestUrl = '{{ request.build_absolute_uri }}';

  for(var i = 0; i < scriptTags.length; i++) {
    var scriptTag = scriptTags[i];

    // src matches the url of this request, and not processed it yet.
    if (scriptTag.src == thisRequestUrl && processedScripts.indexOf(scriptTag) < 0) {
      processedScripts.push(scriptTag);

      // add the style tag into the head (once only)
      if(styleTags.length == 0) {
        // add a style tag to the head
        var styleTag = document.createElement("link");
        styleTag.rel = "stylesheet";
        styleTag.type = "text/css";
        styleTag.href =  "http://{{request.get_host}}/static/css/widget.css";
        document.getElementsByTagName('head')[0].appendChild(styleTag);
        styleTags.push(styleTag);
      }

      // Create a div
      var div = document.createElement('div');
      div.id = 'fpbase-{{object.slug}}-widget';

      // add the cleanslate classs for extreme-CSS reset.
      div.className = 'fpbase widget';

      scriptTag.parentNode.insertBefore(div, scriptTag);

      div.innerHTML = '<img width="100%" src="http://{{request.get_host}}/spectra_img/{{object.slug}}.svg?title=1">';

    }
  }
})(this);

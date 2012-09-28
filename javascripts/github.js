var github = (function(){
  function splitTag(tag){
    pos = tag.lastIndexOf("v");
    tag = tag.slice(pos+1);
    taga = tag.split(".");
    return {major: parseInt(taga[0]), minor: parseInt(taga[1]), revision: parseInt(taga[2])}
  }
  function getTag(ref){
    pos = ref.lastIndexOf("/v");
    return ref.slice(pos+1);
  }
  function renderList(options, tags){
    var i = 0, fragment = '', t = $(options.target)[0];

    for(i = 0; i < tags.length; i++) {
      fragment += '<li class="loading"><span style="font-weight:bolder;">'+tags[i].tag+'</span>&nbsp;<a href="https://github.com/'+options.user+'/'+options.repository+'/zipball/'+tags[i].tag+'">zip</a>&nbsp;<a href="https://github.com/'+options.user+'/'+options.repository+'/tarball/'+tags[i].tag+'">tar.gz</a><p id="gh_tag_'+tags[i].sha+'">Additional information not available.</p></li>';
    }
    t.innerHTML = fragment;
    for(i = 0; i < tags.length; i++) {
      var tag = tags[i];
      $.ajax({
          url: "https://api.github.com/repos/"+options.user+"/"+options.repository+"/git/tags/"+tags[i].sha+"?callback=?"
        , type: 'jsonp'
        , success: function(data) {
          if (!data || !data.data || !data.data.sha) return;
          var t = $('#gh_tag_'+data.data.sha)[0];
          var fragment = '';
          if (data.data.tagger && data.data.tagger.date) {
            fragment += '<span style="font-style:oblique;">Date:</span>&nbsp;' + new Date(data.data.tagger.date).toLocaleDateString()
          }
          if (data.data.tagger && data.data.tagger.date && data.data.message) {
            fragment += '<br />';
          }
          if (data.data.message) {
            fragment += '<span style="font-size:0.875em;">' + data.data.message + '</span>';
          }
          if (fragment != '') {
            t.innerHTML = fragment;
          }
	}
      });
    }
  }
  return {
    showTags: function(options){
      $.ajax({
          url: "https://api.github.com/repos/"+options.user+"/"+options.repository+"/git/refs/tags?callback=?"
        , type: 'jsonp'
        , error: function (err) { $(options.target + ' li.loading').addClass('error').text("Error loading feed"); }
        , success: function(data) {
          var tags = [];
          if (!data || !data.data ) { return; }
          for (var i = 0; i < data.data.length; i++) {
            var d = data.data[i];
            if (d.ref.substr(0, 11) == "refs/tags/v") {
              var tag = getTag(d.ref);
              tags.push({ref: d.ref, tag: tag, url: d.object.url, sha: d.object.sha, version: splitTag(tag)});
            }
          }
          tags.sort(function(a, b) {
            var av = a.version,
                bv = b.version;
            if (av.major < bv.major) { return -1; }
            if (av.major > bv.major) { return  1; }
            // av.major == bv.major
            if (av.minor < bv.minor) { return -1; }
            if (av.minor > bv.minor) { return  1; }
            // av.major == bv.major  && av.minor == bv.minor
            if (av.revision < bv.revision) { return -1; }
            if (av.revision > bv.revision) { return  1; }
            return 0;
          });
          tags.reverse();

          if (options.count) { tags.splice(options.count); }
          else { tags.splice(10); } // Do not show more than 10 tags
          renderList(options, tags);
        }
      });
    }
  };
})();

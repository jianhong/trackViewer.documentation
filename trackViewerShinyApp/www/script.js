// Make sure the Shiny connection is established
$(document).on('shiny:connected', function(event) {

  /********** NON-REACTIVE DOM MANIPULATION **********/
  // Detect input change and change UI to reflect that
  $(document).on('shiny:inputchanged', function(event) {
    if (/^file\d+$/.exec(event.name)) {
      fileIndex = Number(event.name.replace("file", ""));
      ext = event.value.substr( (event.value.lastIndexOf('.') +1) );
      var $select = $('#format'+fileIndex).selectize();
      var selectize = $select[0].selectize;
      switch(ext.toLowerCase()){
        case 'bed':
          selectize.setValue(selectize.search("BED").items[0].id);
            break;
        case 'bg':
        case 'bedgraph':
          selectize.setValue(selectize.search("bedGraph").items[0].id);
            break;
        case 'bw':
        case 'bigwig':
            selectize.setValue(selectize.search("BigWig").items[0].id);
            break;
      }
      //$('#sample'+fileIndex).val(event.value.replace("."+ext, ""));
    }
    if (/^lollifile\d+$/.exec(event.name)) {
      fileIndex = Number(event.name.replace("lollifile", ""));
      ext = event.value.substr( (event.value.lastIndexOf('.') +1) );
      var $select = $('#lolliformat'+fileIndex).selectize();
      var selectize = $select[0].selectize;
      switch(ext.toLowerCase()){
        case 'bed':
          selectize.setValue(selectize.search("BED").items[0].id);
            break;
        case 'bg':
        case 'bedgraph':
          selectize.setValue(selectize.search("bedGraph").items[0].id);
            break;
        case 'gz':
        case 'vcf':
            selectize.setValue(selectize.search("VCF").items[0].id);
            break;
        case 'csv':
            selectize.setValue(selectize.search("pie.stack.csv").items[0].id);
            var $type = $('#lollitype'+fileIndex).selectize();
            var typesel = $type[0].selectize;
            typesel.setValue(typesel.search("pie.stack").items[0].id);
            break;
      }
      //$('#lollisample'+fileIndex).val(event.value.replace("."+ext, ""));
    }
  });
  
});
function toggleIcon(e) {
    //console.log(e);
    elem = $(e.target)
        .prev('.panel-heading')
        .find(".more-less");
    //console.log(elem);
    elem.toggleClass('glyphicon-plus glyphicon-minus');
}

function addtableClass() {
    $("table:not(.table)").addClass('table');
}

$(document).ready( function() {
    //$("#includedContent").load("/www/test.html");
    
    $('.panel-collapse').on('hidden.bs.collapse', toggleIcon);
    $('.panel-collapse').on('shown.bs.collapse', toggleIcon);
    
    $(".file-format-content").load('file_formats.html')
    
    addtableClass()
});

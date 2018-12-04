var debug = true
var testing = true

function addResultsLinks(msg) {
    var links = document.getElementsByClassName("results_button");
    if (debug) {
        console.log(links);
    }
    
    var i;
    for (i = 0; i < links.length; i++) { 
        links[i].onclick = function() {
            if (debug) {
                console.log(this);
            }
            // send message to shiny
            var message = { results_source:this.id, rand: Math.random() }
            Shiny.onInputChange("js_results_source", message);
            return false;
        }
    }
}

$(document).ready( function() {
    addResultsLinks('initialise')
    
});

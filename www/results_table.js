var debug = true
var testing = true

//function setParameters(user_data) {
//    if (debug) {
//        console.log(debug, testing);
//        console.log(user_data);
//    }
//    debug = user_data.debug;
//    testing = user_data.testing;
//    
//    if (debug) {
//        console.log(debug, testing);
//    }
//}

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
    
    //Shiny.addCustomMessageHandler('set_up', setParameters)
    //
    //Shiny.addCustomMessageHandler('counts_link', addCountPlotLinks)
});

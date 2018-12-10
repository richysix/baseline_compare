var debug = true
var testing = true

function hightlightSelectedButton(buttonId) {
    var buttons = document.getElementsByClassName("results-button");
    var i;
    for (i = 0; i < buttons.length; i++) {
        if (debug) {
            console.log(buttons[i])
        }
        id = buttons[i].id
        classNames = buttons[i].className.split(" ");
        if (classNames.indexOf('selected-results') == -1) {
            if (buttonId == id) {
                buttons[i].className += " " + 'selected-results';
            }
        } else {
            if (buttonId != id) {
                buttons[i].className = buttons[i].className.replace(/\bselected-results\b/g, "");
            }
        }
        if (debug) {
            console.log(buttons[i])
        }
    }
}

$(document).ready( function() {
    
    Shiny.addCustomMessageHandler("selected_results_button",
                                  hightlightSelectedButton)
});

var debug = true
var testing = true
var sigLevel = 0.05

// function to hightlight which results button has been clicked on
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

// function to update the sigLevel variable based on selection in the
// Shiny UI
function updateSigLevel(msg) {
    if (debug) {
        console.log(msg)
        console.log(sigLevel)
    }
    sigLevel = parseFloat(msg);
    if (debug) {
        console.log(sigLevel)
    }
}

//function to change the colour of table cells base on their content
//used when producing the results DataTable
var table_def = {
    'expt_only': { 'padj': 3, 'log2fc': 2 },
    'plus_baseline': { 'padj': 5, 'log2fc': 4 },
    'with_stage': { 'padj': 7, 'log2fc': 6 },
}

function colourCells( row, data, dataIndex ) {
    for (x in table_def){
        p_idx = table_def[x]['padj']
        l_idx = table_def[x]['log2fc']
        if (data[p_idx] == 'NA'){
            $('td:eq(' + p_idx + ')', row).addClass( 'table-cell-na' );
            $('td:eq(' + l_idx + ')', row).addClass( 'table-cell-na' );
        } else if (data[p_idx] >= sigLevel) {
            $('td:eq(' + p_idx + ')', row).addClass( 'table-cell-notSig' );
            $('td:eq(' + l_idx + ')', row).addClass( 'table-cell-notSig' );
        }
    }
}

$(document).ready( function() {
    
    Shiny.addCustomMessageHandler("selected_results_button",
                                  hightlightSelectedButton)

    Shiny.addCustomMessageHandler("sigLevel",
                                  updateSigLevel)
});

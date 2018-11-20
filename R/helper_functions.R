library(scales)

# set default colour palette
colour_blind_palette <- 
  c( 'blue' = rgb(0,0.45,0.7),
     'yellow' = rgb(0.95, 0.9, 0.25),
     'vermillion' = rgb(0.8, 0.4, 0),
     'purple' = rgb(0.8, 0.6, 0.7),
     'blue_green' = rgb(0, 0.6, 0.5),
     'sky_blue' = rgb(0.35, 0.7, 0.9),
     'black' = rgb(0, 0, 0),
     'orange' = rgb(0.9, 0.6, 0),
     'grey' = 'grey30'
  )

#' Create a colour palette
#'
#' \code{colour_palette} takes a DESeq2DataSet object and creates a colour palette for the stage attribute
#'
#'    If the number of levels of the stage variable is more than available in the
#'    colour blind palette, the hue_pal function from the scales library is used
#'    
#' @param dds        DESeq2DataSet - colData must contain a column named 'stage'
#' 
#' @return A named vector of colours for the levels of stage
#'
#' @examples
#' colour_palette(dds)
#' 
colour_palette <- function(dds) {
  # check number of levels
  num_colours <- nlevels(colData(dds)[['stage']])
  if (num_colours > length(colour_blind_palette)) {
    ord1 <- seq(1,num_colours,2)
    ord2 <- seq(2,num_colours,2)
    colour_palette <- hue_pal()(num_colours)[ order(c(ord1,ord2)) ]
    names(colour_palette) <- levels(colData(dds)[['stage']])
  } else {
    colour_palette <- colour_blind_palette[seq_len(num_colours)]
    names(colour_palette) <- levels(colData(dds)[['stage']])
  }
  return(colour_palette)
}

#' Create a shape palette
#'
#' \code{shape_palette} takes a DESeq2DataSet object and creates a palette of shapes
#'
#'    If the number of levels of the condition variable is more than shapes 
#'    available (5) then an error is thrown. The baseline level is always a circle.
#'    
#' @param dds DESeq2DataSet - colData must contain a column named 'condition'
#' 
#' @return A named vector of shape indices for the levels of condition
#'
#' @examples
#' shape_palette(dds)
#' 
shape_palette <- function(dds) {
  # check number of levels
  num_shapes <- nlevels(colData(dds)[['condition']])
  shapes <- 21:25
  if (num_shapes > length(shapes) + 1) {
    # error message
    stop('There are two many levels of the condition variable to show with shapes')
  } else {
    shape_palette <- shapes[ seq_len(num_shapes) ]
    names(shape_palette) = c('baseline', setdiff(levels(colData(dds)[['condition']]), c('baseline')))
    return(shape_palette)
  }
}



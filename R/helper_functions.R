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
#' \code{colour_palette} takes a factor and creates a colour palette
#'
#'    If the number of levels of the factor is more than available in the
#'    colour blind palette, the hue_pal function from the scales library is used
#'    
#' @param x factor
#' 
#' @return A named vector of colours for the levels of x
#'
#' @examples
#' colour_palette(colData(dds)$stage)
#' 
colour_palette <- function(x) {
  # check this is a factor
  if (class(x) != 'factor') {
    stop('Not a factor!')
  }
  # check number of levels
  num_colours <- nlevels(x)
  if (num_colours > length(colour_blind_palette)) {
    ord1 <- seq(1,num_colours,2)
    ord2 <- seq(2,num_colours,2)
    colour_palette <- hue_pal()(num_colours)[ order(c(ord1,ord2)) ]
    names(colour_palette) <- levels(x)
  } else {
    colour_palette <- colour_blind_palette[seq_len(num_colours)]
    names(colour_palette) <- levels(x)
  }
  return(colour_palette)
}

#' Create a shape palette
#'
#' \code{shape_palette} takes a factor and creates a palette of shapes
#'
#'    If the number of levels of the factor is more than shapes 
#'    available (5) then an error is thrown. 
#'    Also throws an error is one of the levels is not "baseline"
#'    The "baseline" level is always a circle.
#'    
#' @param x factor
#' 
#' @return A named vector of shape indices for the levels of x
#'
#' @examples
#' shape_palette(colData(dds)$condition)
#' 
shape_palette <- function(x) {
  # check this is a factor
  if (class(x) != 'factor') {
    stop('Not a factor!')
  }
  # check for "baseline" level
  if (!('baseline' %in% x)) {
    stop('None of the levels are "baseline"')
  }
  # check number of levels
  num_shapes <- nlevels(x)
  shapes <- 21:25
  if (num_shapes > length(shapes)) {
    # error message
    stop('There are two many levels of the condition variable to show with shapes')
  } else {
    shape_palette <- shapes[ seq_len(num_shapes) ]
    names(shape_palette) = c('baseline', setdiff(levels(x), c('baseline')))
    return(shape_palette)
  }
}

#' Open a graphics device
#'
#' \code{open_graphics_device} takes a file extension and opens a graphics device
#'
#'    The device name must be one of pdf, eps, svg and png. 
#'    Anything else causes a pdf device to be opened
#'    
#' @param dev_name character - file extension, used to determine device type
#' @param filename character - file name
#' 
#' @return NULL (invisibly)
#'
#' @examples
#' open_graphics_device('pdf', 'plot.pdf')
#' 
#' open_graphics_device('svg', 'plot.svg')
#' 
open_graphics_device <- function(dev_name, filename){
  if (dev_name == "pdf") {
    pdf(file = filename, paper = "special", height = 7, width = 10) # open the pdf device
  } else if (dev_name == "eps") {
    postscript(file = filename, paper = "special", height = 7, width = 10) # open the postscript device
  } else if (dev_name == "svg") {
    svglite(file = filename, height = 7, width = 10) # open the svg device
  } else if (dev_name == "png") {
    png(filename = filename, height = 600, width = 960, res = 100) # open the png device
  } else {
    pdf(file = filename, paper = "special", height = 7, width = 10) # open the pdf device
  }
  invisible(NULL)
}

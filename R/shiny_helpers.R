#' Create an invalid column alert 
#'
#' \code{create_invalid_col_alert} takes a column type and creates an alert
#'
#' @param column_type character - column_type 'condition' or 'sex'
#' @param column_data vector - invalid column
#' @param allowed_values vector - allowed values for column
#' 
#' @return NULL (invisibly)
#'
#' @examples
#' create_invalid_col_alert('condition', c('het', 'hom', 'hemi'),
#'                          c('hom', 'het', 'wt', 'mut', 'sib'))
#' 
#' create_invalid_col_alert('sex', c('M', 'F', 'MALE'), c('F', 'M'))
#' 
create_invalid_col_alert <- function(column_type, column_data, allowed_values,
                                     session){
  alert_id <- paste('invalid', column_type, 'col', sep = "_")
  alert_title <- paste('Invalid', 
                       paste0(toupper(substring(column_type, 1, 1)),  substring(column_type, 2)),
                       'Column')
  closeAlert(session, alert_id)
  createAlert(
    session, anchorId = 'invalid_col_alert', dismiss = FALSE,
    alertId = alert_id, title = alert_title,
    content = paste('The selected', column_type,  'column contains', 
                    'entries that not allowed.', 'Valid values are: ',
                    paste0(allowed_values, collapse = ', '), '<br>',
                    'Selected column looks like this: ',
                    paste0(column_data, collapse = ', ')),
    style = 'danger'
  )
  invisible(NULL)
}

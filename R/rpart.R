## implementation of rpart for competing risks data

#' @examples
#' ulist <- list(eval = efun, split = sfun, init = ifun)
#' fit   <- rpart(formula, data, method = ulist, minsplit = 10)
#' 

# efun <- function(y, offset, parms, wt) {
#   list(label, deviance)
# }
# 
# ifun <- function(y, wt, parms) {
#   list(y, parms = NULL, numresp, numy, summary = NULL, text = NULL)
# }
# 
# sfun <- function(y, wt, x, parms, continuous) {
#   list(goodness, direction)
# }

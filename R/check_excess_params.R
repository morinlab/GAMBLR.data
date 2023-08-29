#' @title Check Excess Params
#'
#' @description Function for checking excessive parameter names.
#' This function will notify the user if any unavailable parameters are called for any given given function.
#' This function is designed to work as internal function-call in already available GAMBLR functions.
#'
#' @details Catch function calls containing unsupported arguments.
#'
#' @param ... Parameters to check.
#'
#' @return Nothing
#'
#' @export
#'
check_excess_params = function(...){
  callingFun = as.list(sys.call(-1))[[1]]
  arguments <- list(...)
  extraneous = names(arguments)
  if(length(arguments)>0){
    stop(paste(" You have given one or more unsupported or deprecated argument to ",callingFun,". Please check the documentation and spelling of your arguments.\nOffending argument(s):", paste(extraneous, collapse = ",")), call. = FALSE)
  }
}
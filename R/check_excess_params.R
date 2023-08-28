
#' Catch function calls containing unsupported arguments
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
check_excess_params = function(...){
  callingFun = as.list(sys.call(-1))[[1]]
  arguments <- list(...)
  extraneous = names(arguments)
  if(length(arguments)>0){
    stop(paste("You have given one or more unsupported or deprecated argument to ",callingFun,". Please check the documentation and spelling of your arguments.\nOffending argument(s):",paste(extraneous,collapse = ",")),call.=FALSE)
  }
}

#' @title Separate a chromosome region into chunks
#' 
#' @description `region_to_chunks` breaks the input string that stores a chromosome 
#' region to create a list with chromosome number and start and end positions as 
#' separated elements.
#'
#' @param region A single string that stores a chromosome region. Any format like 
#' "chr1:100000-200000", "1:100000-200000", "chr1:100'000-200'000" is possible. 
#'
#' @return A list with length 3 and names "chromosome", "start" and "end.
#' @export
#'
#' @examples
#' region_to_chunks(region = "chr1:100000-200000")
#' 
region_to_chunks = function(region){
  region = unname(region)
  region = gsub(",", "", region)
  #format is chr6:37060224-37151701
  split_chunks = unlist(strsplit(region, ":"))
  chromosome = split_chunks[1]
  startend = unlist(strsplit(split_chunks[2], "-"))
  qstart = startend[1]
  qend = startend[2]
  return(list(chromosome = chromosome, start = qstart, end = qend))
}

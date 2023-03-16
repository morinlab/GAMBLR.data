
#' @title Get genes from one or more gene lists
#'
#' @description Retrieve gene names from bundled lymphoma gene lists. 
#'
#' @details Complete lists of genes described as significantly mutated in large
#' lymphoma studies have been curated and provided with this package. 
#'
#' @param Optional vector specifying one or more lymphoma entities e.g. MCL, DLBCL, BL
#' @param curated_only  Specify FALSE to retrieve all genes or leave default for the curated subset
#'
#' @return
#' 
#' @import dplyr 
#' @export
#'
#' @examples
#' all_dlbcl_genes = get_genes(entities="DLBCL",curated_only=FALSE)
#' all_curated_genes = get_genes()

get_genes = function(entities=c("DLBCL","MCL","BL"),curated_only=TRUE,gene_format="symbol"){
    #construct file name using entity
    r_objects = entities %>% 
      tolower() %>% 
      paste0(.,"_genes_latest") #update this to get the version from the config if specified
    
  all_genes = c()
  for (r_object in r_objects){
    if(!exists(r_object)){
      warning(paste("Skipping entity because object",r_object, "doesn't exist"))
    }
    this_df = base::get(r_object)
    if(curated_only){
      this_df = dplyr::filter(this_df,curated==TRUE) 
    }
    if(gene_format=="symbol"){
      these_genes =pull(this_df,Gene)
    }else if(gene_format=="ensembl"){
      #TO DO: the BL and MCL files need to be updated to support this option
      stop("this option is not yet supported")
      these_genes =pull(this_df,ensembl_gene_id)
    }
    all_genes = c(all_genes,these_genes)

  }
  return(unique(all_genes))
}
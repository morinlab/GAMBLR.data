
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


#' @title Get standardized colours for lymphoid cancers
#'
#' @description Retrieve and visualize standardized colours schemes for lymphoid cancers
#'
#' @details Colours hand picked to represent various common entities and clinical variables relevant for lymphoid cancers.
#'
#' @param show_available Set to TRUE to see what options are available. 
#' @param this_category  Optionally supply one of the available categories to see a subset of the options that are available
#' @param this_group  Optionally supply one of the available groups to see the palette just for this group
#' @param drop_alias Set to FALSE to see the redundant colours with their aliases
#'
#' @return
#' 
#' @import dplyr 
#' @export
#'
#' @examples 
#' get_colours(show_available=TRUE)
#' # printout shows that "subgroup" is one option to narrow it down, supply this to the function as this_category:
#' get_colours(show_available=TRUE,this_category = "subgroup")
#' # printout and plot shows several options. Pick the one you want to visualize it in isolation
#' get_colours(show_available=TRUE,this_group = "LymphGen")
#' #if satisfied, get the result for use with ggplot
#' col_vec = get_colours(this_group = "LymphGen",as_named_vector = T)
#' #ggplot(some_data) + scale_fill_manual(values=col_vec) 

get_colours = function(show_available=FALSE,this_category,this_group,as_named_vector=FALSE,drop_alias=TRUE){
  if(drop_alias){
    colour_codes = dplyr::filter(colour_codes,is.na(is_alias))
  }
  if(show_available){
    if(missing(this_category) & missing(this_group)){
      print("Supply a category using this_category parameter. Current options for 'category' are:")
      p = colour_codes  %>% 
        group_by(category) %>% 
        mutate(n=n()) %>% 
        slice_head() %>%
        rename(c("example"="name"))
      print(p)
      print("Supply a group using this_group parameter. Current options for 'group' are:")
      p = colour_codes  %>% 
        group_by(group) %>% 
        mutate(n=n()) %>% 
        slice_head() %>%
        rename(c("example"="name"))
      print(p)
    }else if(!missing(this_group)){
      p = colour_codes %>% 
        dplyr::filter(group == this_group) %>%
        ggplot(aes(x=name,fill=name,y=1)) + geom_col() + 
        theme(legend.position="none") + facet_wrap(~category,scales = "free_y") + 
        scale_fill_manual(values=allcols) + coord_flip()
      print(p)
      return()
    }else if(!missing(this_category)){
      p = colour_codes %>% dplyr::filter(category==this_category) %>% 
        group_by(group) %>% mutate(n=n()) %>% slice_head()
      message("Current options for 'group' within this category are:")
      g = pull(p,group)
      g = paste(g,collapse=",")
      message(g)
      p = colour_codes %>% 
        dplyr::filter(category==this_category) %>%
        ggplot(aes(x=name,fill=name,y=1)) + geom_col() + 
        theme(legend.position="none") + facet_wrap(~group,scales = "free_y") + 
        scale_fill_manual(values=allcols) + coord_flip()
      print(p)
      return()
      
    }
    
  }
  if(missing(this_category) & missing(this_group)){
    stop("provide a category or group via this_category or this_group. To see what's available run this function with show_available = TRUE")
  }else if (!missing(this_group)){
    colour_list = dplyr::filter(colour_codes,group==this_group) %>%
      dplyr::select(name,colour) %>% 
      column_to_rownames("name")
  }else if (!missing(this_category)){
    colour_list = dplyr::filter(colour_codes,category==this_category) %>%
      dplyr::select(name,colour) %>% 
      column_to_rownames("name")
  }
  if(as_named_vector){
    #useful for ggplot scale_X_manual
    allcols = colour_list$colour
    names(allcols) = rownames(colour_list)
    return(allcols)
  }
  return(colour_list)
}

#' @title Get CNV and coding SSM combined status
#' 
#' @description For each specified chromosome region (gene name), return status 1 if the copy number (CN) 
#'   state is non-neutral, *i.e.* different from 2, or if the region contains any coding simple somatic mutation (SSM).
#' 
#' @details The user can choose from which regions are intended to return only copy number variation (CNV) status, 
#'   only coding SSM status, or at least the presence of one of them. This behavior is controlled by the arguments 
#'   `genes_and_cn_threshs` (column `cn_thresh`) and `only_cnv`.
#'   
#'   This function internally calls the `get_cn_states`, `get_ssm_by_samples` and `get_coding_ssm_status`functions. 
#'   Therefore, many of its arguments are assigned to these functions. If needed, see the documentation of these 
#'   functions for more information. 
#'   
#'   In the case of returning NA values, this is due to the `get_cn_segments` function not being able to internally 
#'   return any copy number segments from the specified chromosome region.
#' 
#' @param genes_and_cn_threshs A data frame with columns "gene_id" and "cn_thresh". The "gene_id" column stores 
#'   gene symbols (characters) which determine the regions to return CNV and/or coding SSM status. The "cn_thresh" 
#'   column stores integers that mean the maximum or minimum CN states to return status 1 (contains CNV) for 
#'   its respective gene. If this integer is below 2 (neutral CN state for diploids), it is taken as the maximum 
#'   (gene consider as tumor suppressor); if above 2, it is the minimum (oncogene); if equal to 2, do not consider
#'   CNV to return status.
#' @param these_samples_metadata The metadata for samples of interest to be included in the returned matrix.
#'   Can be created with `get_gambl_metadata` function.
#' @param this_seq_type The seq type to get results for. Possible values are "genome" (default) or "capture".
#' @param only_cnv A vector of gene names indicating the genes for which only CNV status should be considered, 
#'   ignoring SSM status. Set this argument to "all" or "none" (default) to apply this behavior to all or none 
#'   of the genes, respectively.
#' @param genome_build Reference genome build. Possible values are "grch37" (default) or "hg38".
#' @param from_flatfile Logical parameter indicating whether to use flat file to retrieve mutations. Set to FALSE 
#' to use database instead. Default is TRUE.
#' @param include_hotspots Logical parameter indicating whether hotspots object should also be tabulated. Default is TRUE.
#' @param review_hotspots Logical parameter indicating whether hotspots object should be reviewed to include 
#'   functionally relevant mutations or rare lymphoma-related genes. Default is TRUE.
#' @param subset_from_merge Argument to internally pass to `get_ssm_by_samples` function. If set to TRUE, 
#'   the data will be subset from a pre-merged MAF of samples with the specified this_seq_type, Instead of merging 
#'   individual MAFs. Default is FALSE.
#' @param augmented Argument to internally pass to the functions `get_ssm_by_samples` and `get_coding_ssm_status`. 
#'   A logical parameter (default: TRUE). Set to FALSE to use multi-sample patients, instead of the original MAF 
#'   from each sample.
#' @param min_read_support_ssm Only consider SSMs with at least this many reads in t_alt_count (for cleaning 
#'   up augmented MAFs).
#'
#' @return A data frame with CNV and SSM combined status.
#' 
#' @import dplyr
#' @export
#'
#' @examples
#' # Define samples
#' these_sample_ids = c(
#'   "BLGSP-71-06-00160-01A-03D",
#'   "BLGSP-71-06-00252-01A-01D",
#'   "BLGSP-71-19-00122-09A.1-01D",
#'   "BLGSP-71-19-00523-09A-01D",
#'   "BLGSP-71-21-00187-01A-01E",
#'   "BLGSP-71-21-00188-01A-04E"
#' )
#' 
#' # Get sample meta data
#' this_meta = get_gambl_metadata()
#' this_meta = dplyr::filter(this_meta, sample_id %in% these_sample_ids)
#' 
#' # For MYC and SYNCRIP, return CNV and SSM combined status; for MIR17HG, 
#' # return only CNV status; for CCND3 return only SSM status
#' genes_and_cn_threshs = data.frame(
#'   gene_id=c("MYC", "MIR17HG", "CCND3", "SYNCRIP"),
#'   cn_thresh=c(3, 3, 2, 1)
#' )
#' get_cnv_and_ssm_status(
#'   genes_and_cn_threshs,
#'   this_meta,
#'   only_cnv = "MIR17HG",
#' )
#' 
#' # For all genes, return only CNV status
#' genes_and_cn_threshs = data.frame(
#'   gene_id=c("MYC", "MIR17HG", "SYNCRIP"),
#'   cn_thresh=c(3, 3, 1)
#' )
#' get_cnv_and_ssm_status(
#'   genes_and_cn_threshs,
#'   this_meta,
#'   only_cnv = "all",
#' )
#' 
get_cnv_and_ssm_status = function(genes_and_cn_threshs,
                                  these_samples_metadata,
                                  this_seq_type = "genome",
                                  only_cnv = "none",
                                  genome_build = "grch37",
                                  from_flatfile = TRUE,
                                  include_hotspots = TRUE,
                                  review_hotspots = TRUE,
                                  subset_from_merge = FALSE,
                                  augmented = TRUE,
                                  min_read_support_ssm = 3){
  
  # check parameters
  stopifnot('`genes_and_cn_threshs` argument is missing.' = !missing(genes_and_cn_threshs))
  
  stopifnot('`genes_and_cn_threshs` argument must be a data frame with columns "gene_id" (characters) and "cn_thresh" (integers).' = {
    k = class(genes_and_cn_threshs) == "data.frame" &
      all( c("gene_id", "cn_thresh") %in% names(genes_and_cn_threshs) )
    if(k){
      is.character(genes_and_cn_threshs$gene_id) &
        is.numeric(genes_and_cn_threshs$cn_thresh) &
        identical(genes_and_cn_threshs$cn_thresh, trunc(genes_and_cn_threshs$cn_thresh))
    }
  })
  
  stopifnot('`genome_build` argument must be "grch37" or "hg38."' = genome_build %in% c("grch37", "hg38"))
  
  stopifnot('`this_seq_type` argument must be "genome" or "capture."' = this_seq_type %in% c("genome", "capture"))
  
  stopifnot('`only_cnv` argument must be "none", "all", or a subset of `genes_and_cn_threshs$gene_id`' = {
    only_cnv == "none" |
      only_cnv == "all" |
      all(only_cnv %in% genes_and_cn_threshs$gene_id)
  })
  
  # define variables
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter=this_seq_type)
  }else{
    these_samples_metadata = dplyr::filter(these_samples_metadata, seq_type==this_seq_type)
  }
  
  # get gene regions
  my_regions = GAMBLR.utils::gene_to_region(gene_symbol = genes_and_cn_threshs$gene_id,
                                            projection = genome_build,
                                            sort_regions = FALSE)
  
  if(length(my_regions) < nrow(genes_and_cn_threshs)){
    genes_and_cn_threshs = dplyr::filter(genes_and_cn_threshs, gene_id %in% names(my_regions))
  }
  
  ### cnv
  thresh_2 = genes_and_cn_threshs$cn_thresh == 2
  genes_and_cn_threshs_non_neutral = genes_and_cn_threshs[!thresh_2,]
  check_cnv = nrow(genes_and_cn_threshs_non_neutral) > 0
  
  if(check_cnv){
    # get cn states
    cn_matrix = get_cn_states(
      regions_list = my_regions[genes_and_cn_threshs_non_neutral$gene_id],
      region_names = genes_and_cn_threshs_non_neutral$gene_id,
      these_samples_metadata = these_samples_metadata,
      this_seq_type = this_seq_type
    )
    cn_matrix = cn_matrix[these_samples_metadata$sample_id,, drop=FALSE]
    
    # get cnv status
    cnv_status = mapply(function(cnstate, thresh){
      if(thresh < 2){
        cnstate <= thresh
      }else if(thresh == 2){
        cnstate != thresh
      }else if(thresh > 2){
        cnstate >= thresh
      }
    }, cn_matrix, genes_and_cn_threshs_non_neutral$cn_thresh, USE.NAMES = TRUE, SIMPLIFY = FALSE) %>% 
      as.data.frame %>% 
      {. * 1}
    rownames(cnv_status) = rownames(cn_matrix)
    
    # if only CNV statuses are desired, output them
    if(only_cnv == "all"){
      return(cnv_status)
    }
    
    # add cnv status as zero to genes whose cn threshold is 2
    if(any(thresh_2)){
      cnv_status = genes_and_cn_threshs$gene_id[thresh_2] %>% 
        { matrix(0, nrow = nrow(cnv_status), ncol = length(.), dimnames = list(NULL, .)) } %>% 
        cbind(cnv_status)
    }
    cnv_status = dplyr::select(cnv_status, genes_and_cn_threshs$gene_id)
  }
  
  ### ssm
  # genes to get ssm status
  if(only_cnv == "nome"){
    genes_to_check_ssm = genes_and_cn_threshs$gene_id
  }else{
    genes_to_check_ssm = genes_and_cn_threshs$gene_id [ !(genes_and_cn_threshs$gene_id %in% only_cnv) ]
  }
  
  # get maf data
  my_maf = get_ssm_by_samples(
    these_samples_metadata = these_samples_metadata,
    projection = genome_build,
    this_seq_type = this_seq_type,
    min_read_support = min_read_support_ssm,
    these_genes = genes_to_check_ssm,
    subset_from_merge = subset_from_merge,
    augmented = augmented
  )
  
  # get ssm status
  ssm_status = get_coding_ssm_status(
    gene_symbols = genes_to_check_ssm,
    these_samples_metadata = these_samples_metadata,
    maf_data = my_maf,
    projection = genome_build,
    genome_build = genome_build,
    min_read_support = min_read_support_ssm,
    from_flatfile = from_flatfile,
    include_hotspots = include_hotspots,
    include_silent = FALSE,
    augmented = augmented
  ) %>% column_to_rownames("sample_id")
  
  # add missing regions to ssm_status as zero statuses
  missing_regions = !(genes_and_cn_threshs$gene_id %in% names(ssm_status))
  if(any(missing_regions)){
    ssm_status = genes_and_cn_threshs$gene_id[missing_regions] %>% 
      { matrix(0, nrow(ssm_status), length(.), dimnames = list(NULL, .)) } %>% 
      cbind(ssm_status) 
  }
  ssm_status = dplyr::select(ssm_status, genes_and_cn_threshs$gene_id)
  ssm_status = ssm_status[these_samples_metadata$sample_id,, drop=FALSE]
  
  # combine cnv and ssm status
  if(check_cnv){
    (cnv_status | ssm_status) * 1
  }else{
    ssm_status
  }
}

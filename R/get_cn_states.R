#' @title Get CN States.
#'
#' @description Get a copy number matrix for all samples based on segmented data in the database.
#'
#' @details This function returns CN states for the specified regions.
#' For how to specify regions, refer to the parameter descriptions and function examples.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::assign_cn_to_ssm], [GAMBLR.results::get_cn_segments], [GAMBLR.results::get_sample_cn_segments]
#'
#' @param regions_list A vector of regions in the format chrom:start-end.
#' @param regions_bed A bed file with one row for each region you want to determine the CN state from.
#' @param region_names Subset CN states on specific regions (gene symbols e.g FCGR2B).
#' @param these_samples_metadata A metadata table to auto-subset the data to samples in that table before returning.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param all_cytobands Include all cytobands, default is set to FALSE. Currently only supports hg19.
#' @param use_cytoband_name Use cytoband names instead of region names, e.g p36.33.
#' @param missing_data_as_diploid Fill in any sample/region combinations with missing data as diploid (e.g., CN state like 2). Default is FALSE.
#'
#' @return Copy number matrix.
#'
#' @import dplyr circlize tibble stringr tidyr
#' @export
#'
#' @examples
#' #basic usage, generic lymphoma gene list
#' cn_matrix = get_cn_states(regions_bed=GAMBLR.data::grch37_lymphoma_genes_bed)
#'
#' myc_region <- GAMBLR.utils::gene_to_region(
#'  gene_symbol = "MYC",
#'  projection = "grch37",
#'  return_as = "region"
#' )
#'
#' single_gene_cn <- get_cn_states(
#'  regions_list = myc_region,
#'  region_names = "MYC"
#' )
#'
#' # For capture
#' single_gene_cn <- get_cn_states(
#'  regions_list = myc_region,
#'  region_names = "MYC",
#'  this_seq_type = "capture"
#' )
#'
get_cn_states = function(regions_list,
                         regions_bed,
                         region_names,
                         these_samples_metadata,
                         this_seq_type = "genome",
                         all_cytobands = FALSE,
                         use_cytoband_name = FALSE,
                         missing_data_as_diploid = FALSE){

  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter=this_seq_type)
  }else{
    these_samples_metadata = dplyr::filter(these_samples_metadata,seq_type==this_seq_type)
  }
  if(all_cytobands){
    message("Currently, only grch37 is supported")
  }
  #retrieve the CN value for this region for every segment that overlaps it
  bed2region=function(x){
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }
  if(all_cytobands){
    message("Cytobands are in respect to hg19. This will take awhile but it does work, trust me!")
    use_cytoband_name = TRUE
    regions_bed = circlize::read.cytoband(species = "hg19")$df
    colnames(regions_bed) = c("chromosome_name", "start_position", "end_position", "name", "dunno")
    if(use_cytoband_name){
      regions_bed = mutate(regions_bed, region_name = paste0(str_remove(chromosome_name, pattern = "chr"), name))
      region_names = pull(regions_bed, region_name)
    }else{
      region_names = pull(regions_bed, region_name)
    }
    regions = apply(regions_bed, 1, bed2region)
    #use the cytobands from the circlize package (currently hg19 but can extend to hg38 once GAMBLR handles it) Has this been updated?
  }else if(missing(regions_list)){
    if(!missing(regions_bed)){
      regions = apply(regions_bed, 1, bed2region)
    }else{
      warning("You must supply either regions_list or regions_df")
    }
  }else{
    regions = regions_list
  }
  region_segs = lapply(regions,function(x){get_cn_segments(region = x, streamlined = TRUE, this_seq_type = this_seq_type)})
  if(missing(region_names) & !use_cytoband_name){
    region_names = regions
  }
  tibbled_data = tibble(region_segs, region_name = region_names)
  unnested_df = tibbled_data %>%
    unnest_longer(region_segs)

  seg_df = data.frame(ID = unnested_df$region_segs$ID, CN = unnested_df$region_segs$CN,region_name = unnested_df$region_name)
  #arbitrarily take the first segment for each region/ID combination
  seg_df = seg_df %>%
    dplyr::group_by(ID, region_name) %>%
    dplyr::slice(1) %>%
    dplyr::rename("sample_id" = "ID")

  meta_arranged = these_samples_metadata %>%
    dplyr::select(sample_id, pathology, lymphgen) %>%
    arrange(pathology, lymphgen)

  eg = expand_grid(sample_id = pull(meta_arranged, sample_id), region_name = as.character(unique(seg_df$region_name)))
  all_cn = left_join(eg, seg_df, by = c("sample_id" = "sample_id", "region_name" = "region_name"))

  #fill in any sample/region combinations with missing data as diploid
  if(missing_data_as_diploid){
    all_cn = mutate(all_cn, CN = replace_na(CN, 2))
  }

  cn_matrix = pivot_wider(all_cn, id_cols = "sample_id", names_from = "region_name", values_from = "CN") %>%
    column_to_rownames("sample_id")

  #order the regions the same way the user provided them for convenience
  cn_matrix = cn_matrix[,region_names, drop=FALSE]

  return(cn_matrix)
}

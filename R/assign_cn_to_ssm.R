#' @title Assign CN to SSM.
#'
#' @description Annotate mutations with their copy number information.
#'
#' @details This function takes a sample ID with the `this_sample_id` parameter
#'      and annotates mutations with copy number information. A variety of
#'      parameters are at hand for a customized workflow. For example,
#'      the user can specify if only coding mutations are of interest. To do so,
#'      set `coding_only = TRUE`. This function internally calls
#'      `get_ssm_by_samples` and `get_sample_cn_segments`. This function can
#'      also take a vector with genes of interest (`genes`) that the returned
#'      data frame will be restricted to.
#'
#' @param this_sample_id Sample ID of the sample you want to annotate.
#' @param genes A vector of characters with gene symbols (Hugo).
#' @param this_seq_type Specified seq type for returned data. Default is genome.
#' @param projection Specified genome projection that returned data is in
#'      reference to. Default is grch37.
#' @param coding_only Optional. Set to TRUE to restrict to only coding variants
#'      (ssm). Deafult is FALSE.
#' @param assume_diploid Optional, this parameter annotates every mutation as
#'      copy neutral. Default is FALSE.
#' @param include_silent Logical parameter indicating whether to include silent
#'      mutations into coding mutations. Default is FALSE. This parameter only
#'      makes sense if `coding_only` is set to TRUE.
#' @param ... Any additional parameters.
#'
#' @return A list containing a data frame (MAF-like format) with three extra
#'      columns:
#'      - log.ratio is the log ratio from the seg file (NA when no overlap).
#'      - LOH
#'      - CN (the rounded absolute copy number estimate of the region based on
#'          log.ratio, NA when no overlap was found).
#'
#' @import dplyr
#' @export
#'
#' @examples
#' cn_list = assign_cn_to_ssm(
#'      this_sample_id = "DOHH-2",
#'      coding_only = TRUE
#' )
#'
assign_cn_to_ssm = function(
    this_sample_id,
    genes,
    this_seq_type = "genome",
    projection = "grch37",
    coding_only = FALSE,
    assume_diploid = FALSE,
    include_silent = FALSE,
    ...
){

    #warn/notify the user what version of this function they are using
    message("Using the bundled CN segments (.seg) calls in GAMBLR.data...")

    #check if any invalid parameters are provided
    check_excess_params(...)

    #ensure only one sample ID is provided
    if(length(this_sample_id) > 1){
        stop(
            "This function only supports queries of 1 sample ID at the time..."
        )
    }

    #get maf
    maf_sample = get_ssm_by_sample(
        this_sample_id = this_sample_id,
        projection = projection,
        this_seq_type = this_seq_type
    )

    #maf filtering
    #silent mutations
    if(!include_silent){
        coding_class = coding_class[coding_class != "Silent"]
    }

    #coding mutations
    if(coding_only){
        maf_sample = dplyr::filter(
            maf_sample,
            Variant_Classification %in% coding_class
        )
    }

    #subset to genes of interest
    if(!missing(genes)){
        maf_sample = dplyr::filter(maf_sample, Hugo_Symbol %in% genes)
        if(nrow(maf_sample) == 0){
            stop("No variants left after filtering on the provided genes...")
        }
    }

    #get seg
    seg_sample = get_sample_cn_segments(
        these_sample_ids = this_sample_id,
        projection = projection,
        this_seq_type = this_seq_type
    )

    #annotate all CN segments as copy number neutral
    if(assume_diploid){
        diploid = dplyr::mutate(maf_sample, CN = 2)
        return(list(maf = diploid))
    }

    #wrangle the seg file
    seg_sample = seg_sample %>%
        dplyr::filter(end - start > 100) %>%
        mutate(chrom = gsub("chr", "", chrom)) %>%
        rename(
            Chromosome = chrom,
            Start_Position = start,
            End_Position = end,
            LOH = LOH_flag
        ) %>%
        mutate(across(LOH, as.factor))

    #perform an overlap join and add CN columns from the seg file and subset
    # MAF to basic columns (first 45)
    maf_tmp = cool_overlaps(maf_sample, seg_sample, type = "any")

    #rename and change order of columns to match expected format
    maf_with_segs = maf_tmp %>%
        rename(
            Start_Position = Start_Position.x,
            End_Position = End_Position.x
        ) %>%
        dplyr::select(
            colnames(maf_sample),
            LOH, log.ratio, CN
        )

    return(
        list(
            maf = maf_with_segs,
            seg = seg_sample
        )
    )
}

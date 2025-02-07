#' @title Assign CN to SSM.
#'
#' @description Annotate mutations with their copy number information.
#'
#' @details This function takes a metadata table and returns all mutations
#'      for the samples in that metadata. Each mutation is annotated with the 
#'      local copy number state of each mutated site. The user can specify if
#'      only coding mutations are of interest. To do so,
#'      set `coding_only = TRUE`. When necessary, this function relies on
#'      `get_ssm_by_samples` and `get_cn_segments` to obtain the required data. 
#' @param these_samples_metadata Metadata table with one or more rows to specify
#'      the samples to process.
#' @param maf_data A data frame of mutations in MAF format or maf_data object 
#'      (e.g. from `get_coding_ssm` or `get_ssm_by_sample`).
#' @param seg_data A data frame of segmented copy number data or seg_data object
#' @param projection Specified genome projection that returned data is relative to. 
#'      This is only required when it cannot be inferred from maf_df or seg_df 
#'      (or they are not provided). 
#' @param coding_only Optional. Set to TRUE to restrict to only variants in coding space
#'      Default is to work with genome-wide variants.
#' @param assume_diploid Optional, this parameter annotates every mutation as
#'      copy neutral. Default is FALSE.
#' @param include_silent Logical parameter indicating whether to include silent
#'      mutations in coding space. Default is FALSE. This parameter only
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
#' # long-handed way
#' # 1. get some metadata for a collection of samples
#' some_meta = get_gambl_metadata() %>%
#'         dplyr::filter(cohort=="FL_Dreval",
#'         grepl("SP",sample_id))
#' # 2. Get the SSMs for these samples
#' 
#' ssm_genomes_grch37 = get_coding_ssm(projection = "grch37",
#'                                   these_samples_metadata = some_meta)
#' # peek at the results
#' ssm_genomes_grch37 %>% dplyr::select(1:8)
#' 
#' # 3. Lazily let this function obtain the corresponding seg_data for the right genome_build
#' cn_list = assign_cn_to_ssm(some_meta,ssm_genomes_grch37)
#' 
#' cn_list$maf %>% dplyr::select(1:8,log.ratio,CN)
#' 
#' # This won't work because the hg38 seg_data is not bundled
#' ssm_genomes_hg38 = get_coding_ssm(projection = "hg38",
#'                                   these_samples_metadata = some_meta)
#' cn_list = assign_cn_to_ssm(some_meta,ssm_genomes_hg38)
#'
#' # Easiest/laziest way:
#' cn_list = assign_cn_to_ssm(projection = "grch37")
#' 
#' 
#' cn_list$maf %>% dplyr::group_by(Tumor_Sample_Barcode,CN) %>%
#'   dplyr::count()
#'
assign_cn_to_ssm = function(
    these_samples_metadata,
    maf_data,
    seg_data,
    projection,
    coding_only = FALSE,
    assume_diploid = FALSE,
    include_silent = FALSE,
    ...
){
    if(missing(these_samples_metadata)){
        stop("No metadata provided. these_samples_metadata is required")
    }
    #check if any invalid parameters are provided
    check_excess_params(...)
    genomic_data = list()
    if(!missing(maf_data)){
      genomic_data[["maf_data"]] = maf_data
    }
    if(!missing(seg_data)){
      genomic_data[["seg_data"]] = seg_data
    }

    projection <- check_get_projection(genomic_data, suggested = projection)

    if(missing(seg_data)){
      seg_sample = get_cn_segments(
        these_samples_metadata =  these_samples_metadata,
        projection = projection
      )
      missing_from_seg = dplyr::filter(these_samples_metadata,
                                  !sample_id %in% seg_sample$ID) %>% 
        pull(sample_id) %>%
        unique()
      if(length(missing_from_seg) == length(unique(these_samples_metadata$sample_id))){
        stop(paste("No seg_data could be found for ANY of the samples provided for",projection))
      }
      if(length(missing_from_seg)){
        warning(paste("missing seg_data for",length(missing_from_seg),"samples"))
      }
    }else{
      seg_sample = seg_data
    }
    
    if(missing(maf_data)){
      #get maf
      maf_sample = get_ssm_by_samples(
        these_samples_metadata = these_samples_metadata,
        projection = projection,
      )
      missing_from_maf = dplyr::filter(these_samples_metadata,
                                       !sample_id %in% maf_sample$Tumor_Sample_Barcode) %>% 
        pull(sample_id) %>%
        unique()
      if(length(missing_from_maf) == length(unique(these_samples_metadata$sample_id))){
        stop(paste("No mutation could be found for ANY of the samples provided for",projection))
      }
      if(length(missing_from_maf)){
        warning(paste("missing mutation for",length(missing_from_maf),"samples"))
      }
    }else{
      maf_sample = maf_data
    }

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



    #annotate all CN segments as copy number neutral
    if(assume_diploid){
        diploid = dplyr::mutate(maf_sample, CN = 2)
        return(list(maf = diploid))
    }

    #wrangle the seg file
    seg_sample = seg_sample %>%
        dplyr::filter(end - start > 100) %>%
        rename(
            Chromosome = chrom,
            Start_Position = start,
            End_Position = end,
            LOH = LOH_flag,
            Tumor_Sample_Barcode = ID
        ) %>%
        mutate(across(LOH, as.factor))
   
    #perform an overlap join and add CN columns from the seg file and subset
    # MAF to basic columns (first 45)
    maf_tmp = cool_overlaps(maf_sample, seg_sample, 
                            type = "any",
                            columns1=c("Chromosome","Start_Position","End_Position","Tumor_Sample_Barcode"),
                            columns2=c("Chromosome","Start_Position","End_Position","Tumor_Sample_Barcode"))

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

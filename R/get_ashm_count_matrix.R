#' @title Get ASHM Count Matrix.
#'
#' @description Prepare a matrix with one row per sample and one column per
#' region using a set of hypermutated regions.
#'
#' @details Values are the number of mutations in that patient in the region.
#'
#' @param regions_bed A bed file with one row for each region.
#' @param these_samples_metadata This is used to complete your matrix. All GAMBL
#'      samples will be used by default. Provide a data frame with at least
#'      sample_id for all samples if you are using non-GAMBL data.
#' @param this_seq_type The seq type to return results for. Only used if no
#'      metadata is provided with these_samples_metadata.
#' @param projection Which genome build to use for the mutations 
#' (must match the coordinate system your regions to avoid a nonsense result)
#'
#' @return matrix
#'
#' @import dplyr tibble
#' @export
#'
#' @examples
#' regions_bed = create_bed_data(GAMBLR.data::grch37_ashm_regions,
#'                               fix_names="concat",
#'                               concat_cols=c("gene","region"),
#'                               sep="-")
#' my_meta = get_gambl_metadata() %>% dplyr::filter(pathology=="DLBCL")
#' matrix <- get_ashm_count_matrix(
#'      regions_bed = regions_bed,
#'      this_seq_type = "genome"
#' )
#'
#' #this example intentionally fails 
#'  matrix <- get_ashm_count_matrix(regions_bed=regions_bed,this_seq_type = "genome",
#'                             these_samples_metadata = my_meta,
#'                             projection = "hg38")
#' # Error in get_ashm_count_matrix(
#' # Your projection argument does not match the genome_build of regions_bed
#' 
#' # format the name column to include the chromosome coordinates instead of the gene
#' regions_bed = create_bed_data(GAMBLR.data::hg38_ashm_regions,
#'                            fix_names="concat",
#'                            concat_cols=c("chr_name","hg38_start","hg38_end"),
#'                            sep="-")
#'                            
#'  matrix_hg38 <- get_ashm_count_matrix(regions_bed=regions_bed,this_seq_type = "genome",
#'                             these_samples_metadata = my_meta,
#'                             projection = "hg38")
#'
get_ashm_count_matrix = function(
        regions_bed,
        these_samples_metadata,
        this_seq_type,
        projection = "grch37"
    ){
    if(missing(this_seq_type)){
        if(missing(these_samples_metadata)){
            stop(
                "Please supply either the this_seq_type or a metadata from which it can be retrieved"
            )
        }
        this_seq_type <- these_samples_metadata %>%
            pull(seq_type) %>%
            unique()
    }

    if(missing(regions_bed)){
        message(
            "Using aSHM regions in grch37 genome_build as regions_bed"
        )
        if(projection=="grch37"){
          regions_bed <- GAMBLR.data::grch37_ashm_regions %>%
            mutate(name = paste(gene, region, sep = "_")) %>%
            create_bed_data(genome_build = projection)
        }else if(projection=="hg38"){
          regions_bed <- GAMBLR.data::hg38_ashm_regions %>%
            mutate(name = paste(gene, region, sep = "_")) %>%
            create_bed_data(genome_build = projection)
        }else{
          stop(paste("unsupported genome build",projection))
        }
        
    }else{
      if("bed_data" %in% class(regions_bed)){
        if(!get_genome_build(regions_bed)==projection){
          stop(paste("Your genome_build argument does not match the genome_build of regions_bed",get_genome_build(regions_bed),genome_build))
        }
      }
    }

    

    if(missing(these_samples_metadata)){
        all_meta <- get_gambl_metadata(
            seq_type_filter=this_seq_type
        ) %>%
        dplyr::select(sample_id)
    }else{
        all_meta <- these_samples_metadata %>%
            dplyr::select(sample_id)
    }
  
    ashm_maf <- get_ssm_by_regions(
      regions_bed = regions_bed,
      streamlined = TRUE,
      these_samples_metadata = these_samples_metadata,
      use_name_column = TRUE,
      projection = projection
    )
    # Not sure why this was necessary. Possibly because it's also a data.table?
    ashm_maf = strip_genomic_classes(ashm_maf)

    ashm_counted <- ashm_maf %>%
      group_by(sample_id, region) %>%
      tally()

    
    #fill out all combinations so we can get the cases with zero mutations
    eg <- expand_grid(
        sample_id = pull(all_meta, sample_id),
        region = unique(ashm_counted$region)
    )
    all_counts <- left_join(eg, ashm_counted) %>%
        mutate(n = replace_na(n, 0)) %>%
        unique() #not sure where the duplicates are coming from but its annoying

    all_counts_wide <- pivot_wider(
        all_counts,
        id_cols = sample_id,
        names_from = region,
        values_from = n
    ) %>%
        column_to_rownames(var = "sample_id")

    return(all_counts_wide)
}

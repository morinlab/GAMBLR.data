#' @title Review Hotspots.
#'
#' @description Annotate MAF-like data frome with a hot_spot column indicating recurrent mutations.
#'
#' @details This function takes an annotated MAF (with [annotate_hotspots]) and updates an existing column, "hot_spot", in the same data frame.
#' Genes for hotspot review are supplied with the `genes_of_interest` parameter.
#' Currently only a few sets of genes are supported, see parameter description for more information and limitations.
#' The desired genome build can be specified with `genome_build` parameter. Should be the same as the incoming MAF.
#'
#' @param annotated_maf A data frame in MAF format that has hotspots annotated using the function annotate_hotspots().
#' @param genes_of_interest A vector of genes for hotspot review. Currently only FOXO1, MYD88, CREBBP, NOTCH1, NOTCH2, CD79B and EZH2 are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF file. The default is grch37 genome build.
#'
#' @return The same data frame (as given to the `annotated_maf` parameter) with the reviewed column "hot_spot".
#'
#' @import dplyr
#' @export
#'
#' @examples
#' hot_ssms = review_hotspots(annotate_hotspots(get_coding_ssm(this_seq_type = "genome")),
#'                            genes_of_interest = c("CREBBP"))
#'
review_hotspots = function(annotated_maf,
                           genes_of_interest = c("FOXO1", "MYD88", "CREBBP", "NOTCH1", "NOTCH2", "CD79B", "EZH2"),
                           genome_build){
  if(missing(genome_build)){
    if("maf_data" %in% class(annotated_maf)){
      genome_build = get_genome_build(annotated_maf)
      #drop our S3 classes because these additional attributes seem to cause some problems when the data is subsequently munged.
      annotated_maf = strip_genomic_classes(annotated_maf)
    }else{
      stop("genome_build is required")
    }
  }

  # define the list of genes currently supported for review
  supported_genes = c("FOXO1", "MYD88", "CREBBP", "NOTCH1", "NOTCH2", "CD79B", "EZH2")

  # check genome build because CREBBP coordinates are hg19-based or hg38-based

  if (genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
    coordinates = hotspot_regions_grch37
  }else if(genome_build %in% c("hg38", "grch38", "GRCh38")){
    coordinates = hotspot_regions_hg38
  }else{
    stop("The genome build specified is not currently supported. Please provide MAF file in one of the following cordinates: hg19, grch37, hs37d5, GRCh37, hg38, grch38, or GRCh38")
  }
  # check that at least one of the currently supported genes are present
  if (length(intersect(supported_genes, genes_of_interest))==0){
      stop(paste0("Currently only ",  paste(supported_genes, collapse=", "), " are supported. Please specify one of these genes."))
  }
  # notify user that there is limited number of genes currently supported
  if (length(setdiff(genes_of_interest, supported_genes))>0){
      message(strwrap(paste0("Currently only ", paste(supported_genes, collapse=", "),
                             " are supported. By default only these genes from the supplied list will be reviewed. Reviewing hotspots for genes ",
                             paste(intersect(supported_genes, genes_of_interest), collapse = ", "), ", it will take a second ...")))
  }
  if("FOXO1" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "FOXO1" &
                                        HGVSp_Short == "p.M1?",
                                        "TRUE", hot_spot))
  }

  if("CREBBP" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CREBBP" &
                                        Start_Position > coordinates["CREBBP", "start"] &
                                        End_Position < coordinates["CREBBP", "end"] &
                                        Variant_Classification == "Missense_Mutation",
                                        "TRUE", hot_spot))
  }
  if("EZH2" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "EZH2" &
                                        Start_Position > coordinates["EZH2", "start"] &
                                        End_Position < coordinates["EZH2", "end"],
                                        "TRUE", hot_spot))
  }
  if("MYD88" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "MYD88" &
                                        HGVSp_Short %in% c("p.L273P", "p.L265P"),
                                        "TRUE", hot_spot))
  }
  if("NOTCH1" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "NOTCH1" &
                                        Start_Position < coordinates["NOTCH1", "start"],
                                        "TRUE", hot_spot))
  }
  if("NOTCH2" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "NOTCH2" &
                                        Start_Position < coordinates["NOTCH2", "start"],
                                        "TRUE", hot_spot))
  }

  if("CD79B" %in% genes_of_interest){
      truncating_variants = c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Region", "Splice_Site")
      annotated_maf = annotated_maf %>%
         dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CD79B" &
                                         Start_Position < coordinates["CD79B_trunc", "start"] &
                                         Variant_Classification %in% truncating_variants,
                                         "TRUE", hot_spot)) %>%
          dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CD79B" &
                                          Start_Position < coordinates["CD79B_NONtrunc", "start"] &
                                          ! Variant_Classification %in% truncating_variants,
                                          "TRUE", hot_spot))
  }
  annotated_maf = create_maf_data(annotated_maf,genome_build)
  
  return(annotated_maf)
}

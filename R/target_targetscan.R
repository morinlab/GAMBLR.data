#' @title Integrate miRNA (MicroRNA) data.
#'
#' @description Add miRNA data to MAF file.
#'
#' @details This function considers the positions of the variants and looks for 
#' mirnas targets that are in the same positions to add their names and sites to
#' them.  
#'
#' @param maf Required. MAF data frame (required columns: Chromosome, Start_Position, End_Position).
#' @param mirna_target Data frame (required columns: Chromosome, Start_Position, End_Position, miRNA, sites)
#' @param projection The genome build projection for the variants you are working with (default is grch37)
#'
#' @return data frame in MAF format with mirna and sites columns.
#'
#' @import dplyr tidyr purrr
#'
#' @examples
#' sample_df = target_targetscan (maf)
#'


target_targetscan <- function(
                              maf,
                              mirna_target,
                              projection = "grch37"
){
    if (projection == "grch37") {
      maf$Chromosome <- gsub("chr", "", maf$Chromosome)
    } else {
      # If there is a mix of prefixed and non-prefixed options
      maf$Chromosome <- gsub("chr", "", maf$Chromosome) 
      maf$Chromosome <- paste0("chr", maf$Chromosome)
    }
    if (missing(mirna_target)){
      mirna_target <- GAMBLR.data::mirna_targetscan
    }
    mirna <- mirna_target %>%
            arrange(
              Chromosome,
              Start_Position,
              End_Position
            )
    maf_mi <- pmap(
                list(
                  maf$Chromosome,
                  maf$Start_Position,
                  maf$End_Position
                ),
                function(
                  chrom,
                  start_pos,
                  end_pos
                ) {
                    match <- mirna %>%
                          filter(
                            (Chromosome == chrom) &
                              (Start_Position <= start_pos &
                                 End_Position >= start_pos) |
                                  (Start_Position > start_pos &
                                     Start_Position <= end_pos
                                   )
                          )
                    if (nrow(match) == 0) {
                      data.frame(
                        Chromosome = chrom,
                        Start_Position = start_pos,
                        End_Position = end_pos,
                        miRNA = NA,
                        sites = NA
                      )
                    } else {
                      match %>%
                        mutate(
                          Chromosome = chrom,
                          Start_Position = start_pos,
                          End_Position = end_pos) %>%
                            select(
                              Chromosome,
                              Start_Position,
                              End_Position,
                              miRNA,
                              sites
                            )
                    }
                }
              )
    mi_match <- bind_rows(maf_mi)
    output <- left_join(
                maf,
                mi_match,
                relationship = "many-to-many"
              ) %>%
                distinct()
    return(output)
}
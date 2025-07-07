#' Oncogenes in grch37 genome build.
#'
#' A data frame with the coordinates of lymphoma oncogenes relative to the grch37 genome build.
#'
#' @format ## `grch37_oncogene`
#' A data frame with 19 rows and 5 columns.
#' \describe{
#'   \item{chrom}{Chromosomes without chr-prefix, 1:22.}
#'   \item{start}{Start coordinate for the specified oncogene.}
#'   \item{end}{End coordinate for the specified oncogene.}
#'   \item{gene}{Lymphoma oncogene.}
#'   \item{entrez}{ENTREZ ID for the specified oncogene.}
#' }
"grch37_oncogene"

#' Oncogenes in hg38 genome build.
#'
#' A data frame with the coordinates of lymphoma oncogenes relative to the hg38 genome build.
#'
#' @format ## `hg38_oncogene`
#' A data frame with 19 rows and 5 columns.
#' \describe{
#'   \item{chrom}{Chromosomes without chr-prefix, 1:22.}
#'   \item{start}{Start coordinate for the specified oncogene.}
#'   \item{end}{End coordinate for the specified oncogene.}
#'   \item{gene}{Lymphoma oncogene.}
#'   \item{entrez}{ENTREZ ID for the specified oncogene.}
#' }
"hg38_oncogene"

#' Chromosome Arms grch37.
#'
#' A data frame with the chromosome arm coordinates in respect to grch37.
#'
#' @format ## `chromosome_arms_grch37`
#' A data frame with 48 rows and 4 columns.
#' \describe{
#'   \item{chromosome}{Chromosomes without chr-prefix, 1:22, X and Y.}
#'   \item{start}{Start coordinates for the specified chromosome arm.}
#'   \item{end}{End coordinates for the specified chromosome arm.}
#'   \item{arm}{Chromosome arm, either p or q.}
#' }
#' @keywords internal
"chromosome_arms_grch37"


#' Chromosome Arms hg38.
#'
#' A data frame with the chromosome arm coordinates in respect to hg38.
#'
#' @format ## `chromosome_arms_hg38`
#' A data frame with 48 rows and 4 columns.
#' \describe{
#'   \item{chromosome}{Chromosomes with chr-prefix, 1:22, X and Y.}
#'   \item{start}{Start coordinates for the specified chromosome arm.}
#'   \item{end}{End coordinates for the specified chromosome arm.}
#'   \item{arm}{Chromosome arm, either p or q.}
#' }
#' @keywords internal
"chromosome_arms_hg38"


#' Double Hit Signature Genes With Weights.
#'
#' A data frame with double hit signature genes (both as ensembl IDs and Hugo symbols) and importance scores.
#'
#' @format ## `dhitsig_genes_with_weights`
#' A data frame with 104 rows and 3 columns.
#' \describe{
#'   \item{Ensembl_ID}{Ensembl IDs as factors, 104 levels.}
#'   \item{ImportanceScore}{Numeric column with importance scores.}
#'   \item{Hugo_Symbol}{Gene symbols in Hugo format as a factor with 104 levels.}
#' }
#' @keywords internal
"dhitsig_genes_with_weights"


#' Genes Blacklist.
#'
#' A tibble with gene symbols (Hugo) that falls within blacklisted regions of the genome.
#'
#' @format ## `gene_blacklist`
#' A tibble with 291 rows.
#' \describe{
#'   \item{Gene}{Genes symbols in Hugo format.}
#' }
#' @keywords internal
"gene_blacklist"


#' grch37 Gene Coordinates.
#'
#' All gene coordinates in respect to grch37.
#'
#' @format ## `grch37_gene_coordinates`
#' A data frame with 63,763 rows and 6 columns.
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene ID}
#'   \item{chromosome}{The chromosome that the gene is residing on}
#'   \item{start}{The start coordinates for the gene}
#'   \item{end}{The end coordinates for the gene}
#'   \item{gene_name}{The gene name}
#'   \item{hugo_symbol}{Gene symbol in Hugo format}
#' }
#' @keywords internal
"grch37_gene_coordinates"


#' Lymphoma Genes (grch37).
#'
#' Lymphoma associated genes in respect to grch37.
#'
#' @format ## `grch37_lymphoma_genes_bed`
#' A data frame with 195 rows and 4 columns.
#' \describe{
#'   \item{chromosome_name}{The chromosome for which the gene is residing on}
#'   \item{start_position}{The start coordinate for the gene}
#'   \item{end_position}{The end coordinate for the gene}
#'   \item{hgnc_symbol}{Gene symbol in Hugo format}
#' }
#' @keywords internal
"grch37_lymphoma_genes_bed"


#' grch37 Partner Genes.
#'
#' Translocation partners for oncogenes in with coordinates in respect to grch37.
#'
#' @format ## `grch37_partners`
#' A data frame with 31 rows and 5 columns.
#' \describe{
#'   \item{chrom}{The chromosome for which the gene is residing on}
#'   \item{start}{The start coordinate for the gene}
#'   \item{end}{The end coordinate for the gene}
#'   \item{gene}{Gene symbol in Hugo format}
#'   \item{entrez}{Entrez ID}
#' }
#' @keywords internal
"grch37_partners"

#' hg38 Gene Coordinates.
#'
#' All gene coordinates in respect to hg38.
#'
#' @format ## `hg38_gene_coordinates`
#' A data frame with 63,763 rows and 6 columns.
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene ID}
#'   \item{chromosome}{The chromosome that the gene is residing on}
#'   \item{start}{The start coordinates for the gene}
#'   \item{end}{The end coordinates for the gene}
#'   \item{gene_name}{The gene name}
#'   \item{hugo_symbol}{Gene symbol in Hugo format}
#' }
#' @keywords internal
"hg38_gene_coordinates"


#' Lymphoma Genes (hg38).
#'
#' Lymphoma associated genes in respect to hg38.
#'
#' @format ## `hg38_lymphoma_genes_bed`
#' A data frame with 195 rows and 4 columns.
#' \describe{
#'   \item{chromosome_name}{The chromosome for which the gene is residing on}
#'   \item{start_position}{The start coordinate for the gene}
#'   \item{end_position}{The end coordinate for the gene}
#'   \item{hgnc_symbol}{Gene symbol in Hugo format}
#' }
#' @keywords internal
"hg38_lymphoma_genes_bed"


#' hg38 Partner Genes.
#'
#' Translocation partners for oncogenes in with coordinates in respect to hg38.
#'
#' @format ## `hg38_partners`
#' A data frame with 31 rows and 5 columns.
#' \describe{
#'   \item{chrom}{The chromsome for which the gene is residing on}
#'   \item{start}{The start coordinate for the gene}
#'   \item{end}{The end coordinate for the gene}
#'   \item{gene}{Gene symbol in Hugo format}
#'   \item{entrez}{Entrez ID}
#' }
#' @keywords internal
"hg38_partners"


#' grch37 Hotspot Regions.
#'
#' Mutation hotspot regions in respect to grch37.
#'
#' @format ## `hotspot_regions_grch37`
#' A data frame with 6 rows and 3 columns.
#' \describe{
#'   \item{chrom}{Chromosome for the described region}
#'   \item{start}{The start coordinate for the region}
#'   \item{end}{The end coordinate for the region}
#' }
#' @keywords internal
"hotspot_regions_grch37"


#' hg38 Hotspot Regions.
#'
#' Mutation hotspot regions in respect to hg38.
#'
#' @format ## `hotspot_regions_hg38`
#' A data frame with 6 rows and 3 columns.
#' \describe{
#'   \item{chrom}{Chromosome for the described region}
#'   \item{start}{The start coordinate for the region}
#'   \item{end}{The end coordinate for the region}
#' }
#' @keywords internal
"hotspot_regions_hg38"


#' Lymphoma Genes Comprehensive.
#'
#' A detailed data frame with lymphoma genes, annotated with evidence from literature and aSHM.
#'
#' @format ## `lymphoma_genes_comprehensive`
#' A data frame with 127 rows and 9 columns.
#' \describe{
#'   \item{ensemble_gene_id}{Gene in Ensemble format}
#'   \item{Gene}{Gene symbol in Hugo format.}
#'   \item{Chapuy}{Boolean flag, TRUE if gene verified by the stated study (Chapuy)}
#'   \item{Reddy}{Boolean flag, TRUE if gene verified by the stated study (Reddy)}
#'   \item{LymphGen}{Boolean flag, TRUE if lymphGen}
#'   \item{curated}{Boolean flag, describing if the gene has been curated or not}
#'   \item{other}{Other support for the described gene}
#'   \item{Lacy}{Boolean flag, TRUE if gene verified by the stated study (Lacy)}
#'   \item{aSHM}{Boolean flag for annotating aSHM}
#' }
#' @keywords internal
"lymphoma_genes_comprehensive"


#' Reddy Genes.
#'
#' Genes identified as significantly mutated in DLBCL by the study of Reddy et al.
#'
#' @format ## `reddy_genes`
#' A data frame with 150 rows and 4 columns.
#' \describe{
#'   \item{Input}{Input}
#'   \item{hgnc_symbol}{The HGNC symbol}
#'   \item{Approved name}{Approved name}
#'   \item{HGNC ID}{HGNC ID}
#' }
#' @keywords internal
"reddy_genes"


#' Target Regions grch37.
#'
#' Target regions in respect to grch37.
#'
#' @format ## `target_regions_grch37`
#' A data frame with 295994 rows and 3 columns.
#' \describe{
#'   \item{chrom}{Chromosome for the descriebd region}
#'   \item{start}{Start coordinate of the region}
#'   \item{end}{End coordiante of the region}
#' }
#' @keywords internal
"target_regions_grch37"


#' Target Regions hg38.
#'
#' Target regions in respect to hg38.
#'
#' @format ## `target_regions_hg38`
#' A data frame with 296453 rows and 3 columns.
#' \describe{
#'   \item{chrom}{Chromosome for the descriebd region}
#'   \item{start}{Start coordinate of the region}
#'   \item{end}{End coordiante of the region}
#' }
#' @keywords internal
"target_regions_hg38"


#' Wright Genes With Weights.
#'
#' Description.
#'
#' @format ## `wright_genes_with_weights`
#' A data frame with 210 rows and 3 columns.
#' \describe{
#'   \item{Ensembl_ID}{Gene in Ensembl ID format}
#'   \item{Hugo_Symbol}{Gene symbol in Hugo format}
#'   \item{Weight_tValue}{Weight Value for the specified gene}
#' }
#' @keywords internal
"wright_genes_with_weights"


#' Default mapping table between mutation type (aka, variant classification) to mutation class
#'
#' A dataset containing the mapping table between genomic mutation type (aka, variant classification) to mutation class.
#' This dataset comes from the g3viz package and was obtained via this URL:
#' https://github.com/morinlab/g3viz/tree/master/data
#'
#' @format A data frame with three columns:
#' \describe{
#'   \item{Mutation_Type}{Mutation type, aka, variant classification}
#'   \item{Mutation_Class}{mutation class}
#'   \item{Short_Name}{short name of mutation type}
#' }
#' @examples
#' mutation.table.df
#' @keywords internal
"mutation.table.df"

#' Mapping table between gene.symbol, uniprot.id, and pfam
#'
#' A dataset containing the mapping table between Hugo symbol, UniProt ID, and
#' Pfam ACC. This dataset comes from the g3viz package and was obtained via this URL:
#' https://github.com/morinlab/g3viz/tree/master/data
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{symbol}{Gene symbol}
#'   \item{uniprot}{UniProt ID}
#'   \item{length}{protein length}
#'   \item{start}{starting position of Pfam domain}
#'   \item{end}{ending position of Pfam domain}
#'   \item{hmm.acc}{Pfam accession number}
#'   \item{hmm.name}{Pfam name}
#'   \item{type}{Pfam type, i.e., domain/family/motif/repeat/disordered/coiled-coil}
#' }
#' @examples
#' hgnc2pfam.df
#' @source Pfam (v31.0) and UniProt
#' @keywords internal
"hgnc2pfam.df"


#' Colour Codes
#'
#' A data frame with colour codes (hex) arranged into different categories, groups.
#'
#' @format ## `colour_codes`
#' A data frame with 221 rows and 5 columns.
#' \describe{
#'   \item{category}{Describes category for any given colour.}
#'   \item{group}{Describes the group for any given colour.}
#'   \item{name}{The name for any given colour.}
#'   \item{colour}{Colour annotated in HEX format.}
#'   \item{is_alias}{Describes if the colour has an alias (yes) or not (NA)}
#' }
#' @keywords internal
"colour_codes"


#' GRCh37 ASHM Regions
#'
#' Aberrant Somatic Hyper Mutation (ASHM) regions in respect to GRCh37.
#'
#' @format ## `grch37_ashm_regions`
#' A data frame with 129 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg19_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg19_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
"grch37_ashm_regions"


#' hg38 ASHM Regions
#'
#' Aberrant Somatic Hyper Mutation (ASHM) regions in respect to hg38.
#'
#' @format ## `hg38_ashm_regions`
#' A data frame with 129 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg38_start}{Start coordinate for region in respect to hg38.}
#'   \item{hg38_end}{End coordinate for region in respect to hg38.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
"hg38_ashm_regions"


#' Lymphoma Genes BL v0.1
#'
#' Genes frequently associated with Burkitt Lymphoma (BL). This is version 0.1.
#'
#' @format ## `lymphoma_genes_bl_v0.1`
#' A data frame withh 128 rows and 11 columns.
#' \describe{
#'   \item{ensembl_gene_id}{Gene ID in ensembl format.}
#'   \item{Gene}{Gene symbol in Hugo format.}
#'   \item{earliest_support_BL}{Pubmeed ID to associated study.}
#'   \item{curated}{Boolean variable annotating if the event is currated or not.}
#'   \item{aSHM_target_DLBCL}{Boolean variable annotating if the event is an aSHM target in DLBCL or not.}
#'   \item{Currator_comments}{Comments from the currator.}
#'   \item{Enrichment_DLBCL_BL}{Annotates the event if it is enriched in BL or DLBCL.}
#'   \item{frequency_BL_Thomas}{The frequency of which the event was reported in the Thomas study.}
#'   \item{frequency_BL_Panea}{The frequency of which the event was reported in the Panea study.}
#'   \item{n_BL_Panea_original}{Total number of mutated tumors as originally reported in Panea study.}
#'   \item{frequency_BL_Panea_original}{Frequency of mutation as originally reported in Panea study.}
#' }
#' @keywords internal
"lymphoma_genes_bl_v0.1"


#' Lymphoma Genes BL v0.2
#'
#' Genes frequently associated with Burkitt Lymphoma (BL). This is version 0.2.
#'
#' @format ## `lymphoma_genes_bl_v0.2`
#' A data frame withh 111 rows and 9 columns.
#' \describe{
#'   \item{Gene}{Gene ID in Hugo format.}
#'   \item{Tier}{Tier of gene in this pathology.}
#'   \item{aSHM}{Whether this gene is associated with aSHM.}
#'   \item{QC}{QC of a particular gene.}
#'   \item{Mean.Variant.Quality}{Mean variant quality.}
#'   \item{citekey}{Alphanumeric representation of the citekey where this gene was first described.}
#'   \item{PMID}{Pubmed ID to associated study.}
#'   \item{MutationEffect}{Annotates the effect of the gene mutation.}
#'   \item{Mutation.PMID}{Pubmed ID to associated study where mutation effect is described.}
#' }
#' @keywords internal
"lymphoma_genes_bl_v0.2"

#' Lymphoma Genes BL Latest
#'
#' Genes frequently associated with Burkitt Lymphoma (BL). This is the most up-to-date version of this dataset.
#'
#' @format ## `lymphoma_genes_bl_v_latest`
#' A data frame withh 111 rows and 9 columns.
#' \describe{
#'   \item{Gene}{Gene ID in Hugo format.}
#'   \item{Tier}{Tier of gene in this pathology.}
#'   \item{aSHM}{Whether this gene is associated with aSHM.}
#'   \item{QC}{QC of a particular gene.}
#'   \item{Mean.Variant.Quality}{Mean variant quality.}
#'   \item{citekey}{Alphanumeric representation of the citekey where this gene was first described.}
#'   \item{PMID}{Pubmed ID to associated study.}
#'   \item{MutationEffect}{Annotates the effect of the gene mutation.}
#'   \item{Mutation.PMID}{Pubmed ID to associated study where mutation effect is described.}
#' }
"lymphoma_genes_bl_v_latest"


#' Lymphoma Genes DLBCL v0.1
#'
#' Genes frequently associated with Diffuse large B cell lymphoma (DLBCL). This is version 0.1.
#'
#' @format ## `lymphoma_genes_dlbcl_v0.1`
#' A data frame with 143 rows and 13 columns.
#' \describe{
#'   \item{ensembl_gene_id}{Gene ID in ensembl format.}
#'   \item{Gene}{Gene symbol in Hugo format.}
#'   \item{Chappuy}{Boolean variable stating if the event is described in the study (Chappuy).}
#'   \item{Reddy}{Boolean variable stating if the event is described in the study (Reddy).}
#'   \item{LymphGen}{Boolean variable stating if the event is a described lymphgen or not.}
#'   \item{curated}{Boolean variable annotating if the event is currated or not.}
#'   \item{other_support}{Variable that annotates the event if there are other support available.}
#'   \item{Lacy}{Boolean variable stating if the event is described in the study (Lacy).}
#'   \item{aSHM}{Boolean varaible annotating if the event is considered an aSHM or not.}
#'   \item{known_hotspots}{Boolean varaible annotating if the event is a known hotspot or not.}
#'   \item{earliest_support}{Pubmeed ID to associated study.}
#'   \item{common_alias}{Variable annotating other common aliases for the event, if such exists.}
#'   \item{noncoding_driver_support}{Boolean variable annotating if the event has noncoding driver support or not.}
#' }
#' @keywords internal
"lymphoma_genes_dlbcl_v0.1"

#' Lymphoma Genes DLBCL v0.2
#'
#' Genes frequently associated with Diffuse large B cell lymphoma (DLBCL). This is version 0.2.
#'
#' @format ## `lymphoma_genes_dlbcl_v0.2`
#' A data frame withh 335 rows and 9 columns.
#' \describe{
#'   \item{Gene}{Gene ID in Hugo format.}
#'   \item{Tier}{Tier of gene in this pathology.}
#'   \item{aSHM}{Whether this gene is associated with aSHM.}
#'   \item{QC}{QC of a particular gene.}
#'   \item{Mean.Variant.Quality}{Mean variant quality.}
#'   \item{citekey}{Alphanumeric representation of the citekey where this gene was first described.}
#'   \item{PMID}{Pubmed ID to associated study.}
#'   \item{MutationEffect}{Annotates the effect of the gene mutation.}
#'   \item{Mutation.PMID}{Pubmed ID to associated study where mutation effect is described.}
#'   \item{MutationEffect.citekey}{Alphanumeric representation of the citekey to associated study where mutation effect is described.}
#'   \item{Mutation.PMID}{Whether this gene is a feature in LymphGen classifier.}
#' }
#' @keywords internal
"lymphoma_genes_dlbcl_v0.2"

#' Lymphoma Genes DLBCL Latest
#'
#' Genes frequently associated with Diffuse large B cell lymphoma (DLBCL).
#' This is the most up-to-date version of this dataset.
#'
#' @format ## `lymphoma_genes_dlbcl_v_latest`
#' A data frame withh 335 rows and 9 columns.
#' \describe{
#'   \item{Gene}{Gene ID in Hugo format.}
#'   \item{Tier}{Tier of gene in this pathology.}
#'   \item{aSHM}{Whether this gene is associated with aSHM.}
#'   \item{QC}{QC of a particular gene.}
#'   \item{Mean.Variant.Quality}{Mean variant quality.}
#'   \item{citekey}{Alphanumeric representation of the citekey where this gene was first described.}
#'   \item{PMID}{Pubmed ID to associated study.}
#'   \item{MutationEffect}{Annotates the effect of the gene mutation.}
#'   \item{Mutation.PMID}{Pubmed ID to associated study where mutation effect is described.}
#'   \item{MutationEffect.citekey}{Alphanumeric representation of the citekey to associated study where mutation effect is described.}
#'   \item{Mutation.PMID}{Whether this gene is a feature in LymphGen classifier.}
#' }
"lymphoma_genes_dlbcl_v_latest"


#' Lymphoma Genes v0.0
#'
#' A comprehenssive resource of genes associated with different types of lymphomas.
#'
#' @aliases lymphoma_genes
#'
#' @format ## `lymphoma_genes_lymphoma_genes_v0.0`
#' A data frame with 196 rows and 12 columns.
#' \describe{
#'   \item{Gene}{Gene symbol in Hugo format.}
#'   \item{DLBCL}{Boolean varaible annotating if the event (genne) is associated with DLBCL.}
#'   \item{FL}{Boolean varaible annotating if the event (genne) is associated with FL.}
#'   \item{BL}{Boolean varaible annotating if the event (genne) is associated with BL.}
#'   \item{MCL}{Boolean varaible annotating if the event (genne) is associated with MCL.}
#'   \item{CLL}{Boolean varaible annotating if the event (genne) is associated with CLL.}
#'   \item{ensembl_geene_id}{Gene ID in ensembl format.}
#'   \item{hgnc_symbol}{Gene symbol in HGNC format.}
#'   \item{LymphGen}{Boolean variable stating if the event is a described lymphgen or not.}
#'   \item{Reddy}{Boolean variable stating if the event is described in the study (Reddy).}
#'   \item{Chappuy}{Boolean variable stating if the event is described in the study (Chappuy).}
#'   \item{entrezgene_id}{Gene ID in entrez fromat.}
#' }
#' @keywords internal
"lymphoma_genes_lymphoma_genes_v0.0"


#' Lymphoma Genes MCL v0.1
#'
#' Genes frequently associated with Mantle cell lymphoma (MCL). This is version 0.1.
#'
#' @format ## `lymphoma_genes_mcl_v0.1`
#' A data frame with 53 rows and 13 columns.
#' \describe{
#'   \item{ensembl_gene_id}{Gene ID in ensembl format.}
#'   \item{Gene}{Gene symbol in Hugo format.v}
#'   \item{curated}{Boolean variable annotating if the event is currated or not.}
#'   \item{Bea}{Boolean variable stating if the event is described in the study (Bea).}
#'   \item{Zhang}{Boolean variable stating if the event is described in the study (Zhang).}
#'   \item{Pararajalingam}{Boolean variable stating if the event is described in the study (Pararajalingam).}
#'   \item{Nadeu}{Boolean variable stating if the event is described in the study (Nadeu).}
#'   \item{Other_support}{Variable that annotates the event if there are other support available.}
#'   \item{appx_overall_freq}{Approximate overall frequency.}
#'   \item{gambl_freq}{Frequency of event described in GAMBL.}
#'   \item{common_alias}{Variable annotating other common aliases for the event, if such exists.}
#'   \item{noncoding_driver_support}{Boolean variable annotating if the event has noncoding driver support or not.}
#'   \item{aSHM}{Boolean varaible annotating if the event is considered an aSHM or not.}
#' }
#' @keywords internal
"lymphoma_genes_mcl_v0.1"

#' Lymphoma Genes MCL v0.2
#'
#' Genes frequently associated with Mantle cell lymphoma (MCL). This is version 0.2.
#'
#' @format ## `lymphoma_genes_mcl_v0.2`
#' A data frame withh 69 rows and 16 columns.
#' \describe{
#'   \item{Gene}{Gene ID in Hugo format.}
#'   \item{curated}{Whether this is a manually curated gene.}
#'   \item{Tier}{Tier of gene in this pathology.}
#'   \item{Bea}{Whether this gene was described in Bea study.}
#'   \item{Zhang}{Whether this gene was described in Zhang study.}
#'   \item{Pararajalingam}{Whether this gene was described in Pararajalingam study.}
#'   \item{Nadeu}{Whether this gene was described in Nadeu study.}
#'   \item{Other_support}{Whether there is a support for this gene other than the studies above.}
#'   \item{appx_overall_freq}{Approximate frequency of mutations at this gene in MCL.}
#'   \item{gambl_freq}{Frequency of mutations at this gene in MCL as inferred in GAMBL.}
#'   \item{common_alias}{Alias name for this gene.}
#'   \item{noncoding_driver_suspect}{Whether there are relevant non-coding mutations at this gene.}
#'   \item{aSHM}{Whether this gene is associated with aSHM.}
#'   \item{QC}{QC of a particular gene.}
#'   \item{CIViC}{Whether this gene is reported in the CIViC knowledge base.}
#'   \item{Earliest_support}{The earlist study to describe this gene to be mutated in MCL.}
#'   \item{citekey}{Alphanumeric representation of the citekey where this gene was first described.}
#' }
#' @keywords internal
"lymphoma_genes_mcl_v0.2"


#' Lymphoma Genes MCL Latest
#'
#' Genes frequently associated with Mantle cell lymphoma (MCL).
#' This is the most up-to-date version of this dataset.
#'
#' @format ## `lymphoma_genes_mcl_v_latest`
#' A data frame withh 69 rows and 16 columns.
#' \describe{
#'   \item{Gene}{Gene ID in Hugo format.}
#'   \item{curated}{Whether this is a manually curated gene.}
#'   \item{Tier}{Tier of gene in this pathology.}
#'   \item{Bea}{Whether this gene was described in Bea study.}
#'   \item{Zhang}{Whether this gene was described in Zhang study.}
#'   \item{Pararajalingam}{Whether this gene was described in Pararajalingam study.}
#'   \item{Nadeu}{Whether this gene was described in Nadeu study.}
#'   \item{Other_support}{Whether there is a support for this gene other than the studies above.}
#'   \item{appx_overall_freq}{Approximate frequency of mutations at this gene in MCL.}
#'   \item{gambl_freq}{Frequency of mutations at this gene in MCL as inferred in GAMBL.}
#'   \item{common_alias}{Alias name for this gene.}
#'   \item{noncoding_driver_suspect}{Whether there are relevant non-coding mutations at this gene.}
#'   \item{aSHM}{Whether this gene is associated with aSHM.}
#'   \item{QC}{QC of a particular gene.}
#'   \item{CIViC}{Whether this gene is reported in the CIViC knowledge base.}
#'   \item{Earliest_support}{The earlist study to describe this gene to be mutated in MCL.}
#'   \item{citekey}{Alphanumeric representation of the citekey where this gene was first described.}
#' }
"lymphoma_genes_mcl_v_latest"


#' Sample Data
#'
#' Sample data bundled as a list of 3 elements. Metadata (data frame) and sample data from two projections (grch37 and hg38),
#' Each projection is organized as a list of 3 elements; maf, seg, and bedpe (all data frames).
#'
#' @format ## `sample_data`
#' A list of three elements.
#' Two elements are lists with three data frames in each list and the last element (metadata) is a data frame.
#' \describe{
#'   \item{meta}{A data frame with metadata.}
#'   \item{grch37}{A list containing 3 data frames; maf, seg, and bedpe. All in respect to grch37.}
#'   \item{hg38}{A list containing 3 data frames; maf, seg, and bedpe. All in respect to hg38.}
#' }
#' @keywords internal
"sample_data"


#' Somatic Hypermutation Locations GRCh37 v0.0
#'
#' A data frame with somatic hypermutation locations in respect to GRCh37, version 0.0.
#'
#' @format ## `somatic_hypermutation_locations_GRCh37_v0.0`
#' A data frame with 76 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg19_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg19_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh37_v0.0"


#' Somatic Hypermutation Locations GRCh37 v0.1
#'
#' A data frame with somatic hypermutation locations in respect to GRCh37, version 0.1.
#'
#' @format ## `somatic_hypermutation_locations_GRCh37_v0.1`
#' A data frame with 78 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg19_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg19_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh37_v0.1"


#' Somatic Hypermutation Locations GRCh37 v0.2
#'
#' A data frame with somatic hypermutation locations in respect to GRCh37, version 0.2.
#'
#' @format ## `somatic_hypermutation_locations_GRCh37_v0.2`
#' A data frame with 88 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg19_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg19_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh37_v0.2"


#' Somatic Hypermutation Locations GRCh37 v0.3
#'
#' A data frame with somatic hypermutation locations in respect to GRCh37, version 0.3.
#'
#' @format ## `somatic_hypermutation_locations_GRCh37_v0.3`
#' A data frame with 101 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg19_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg19_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh37_v0.3"

#' Somatic Hypermutation Locations GRCh37 v0.4
#'
#' A data frame with somatic hypermutation locations in respect to GRCh37, version 0.4.
#'
#' @format ## `somatic_hypermutation_locations_GRCh37_v0.4`
#' A data frame with 129 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg19_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg19_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh37_v0.4"

#' Somatic Hypermutation Locations GRCh37 v0.5
#'
#' A data frame with somatic hypermutation locations in respect to GRCh37, version 0.5.
#'
#' @format ## `somatic_hypermutation_locations_GRCh37_v0.5`
#' A data frame with 130 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg19_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg19_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh37_v0.5"

#' Somatic Hypermutation Locations GRCh37 v0.6
#'
#' A data frame with somatic hypermutation locations in respect to GRCh37, version 0.6.
#'
#' @format ## `somatic_hypermutation_locations_GRCh37_v0.6`
#' A data frame with 127 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg19_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg19_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh37_v0.6"

#' Somatic Hypermutation Locations GRCh37 Latest
#'
#' A data frame with somatic hypermutation locations in respect to GRCh37, the latest version.
#'
#' @format ## `somatic_hypermutation_locations_GRCh37_v_latest`
#' A data frame with 127 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg19_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg19_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh37_v_latest"


#' Somatic Hypermutation Locations GRCh38 v0.0
#'
#' A data frame with somatic hypermutation locations in respect to GRCh38, version 0.0.
#'
#' @format ## `somatic_hypermutation_locations_GRCh38_v0.0`
#' A data frame with 76 rows and 7 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg38_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg38_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#'   \item{name}{Location name.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh38_v0.0"


#' Somatic Hypermutation Locations GRCh38 v0.1
#'
#' A data frame with somatic hypermutation locations in respect to GRCh38, version 0.1.
#'
#' @format ## `somatic_hypermutation_locations_GRCh38_v0.1`
#' A data frame with 77 rows and 7 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg38_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg38_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#'   \item{name}{Location name.}v
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh38_v0.1"


#' Somatic Hypermutation Locations GRCh38 v0.2
#'
#' A data frame with somatic hypermutation locations in respect to GRCh38, version 0.2.
#'
#' @format ## `somatic_hypermutation_locations_GRCh38_v0.2`
#' A data frame with 88 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg38_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg38_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh38_v0.2"


#' Somatic Hypermutation Locations GRCh38 v0.3
#'
#' A data frame with somatic hypermutation locations in respect to GRCh38, version 0.3.
#'
#' @format ## `somatic_hypermutation_locations_GRCh38_v0.3`
#' A data frame with 101 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg38_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg38_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh38_v0.3"


#' Somatic Hypermutation Locations GRCh38 v0.4
#'
#' A data frame with somatic hypermutation locations in respect to GRCh38, version 0.4.
#'
#' @format ## `somatic_hypermutation_locations_GRCh38_v0.4`
#' A data frame with 129 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg38_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg38_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh38_v0.4"

#' Somatic Hypermutation Locations GRCh38 v0.5
#'
#' A data frame with somatic hypermutation locations in respect to GRCh38, version 0.5.
#'
#' @format ## `somatic_hypermutation_locations_GRCh38_v0.5`
#' A data frame with 130 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg38_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg38_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh38_v0.5"

#' Somatic Hypermutation Locations GRCh38 v0.6
#'
#' A data frame with somatic hypermutation locations in respect to GRCh38, version 0.6.
#'
#' @format ## `somatic_hypermutation_locations_GRCh38_v0.6`
#' A data frame with 127 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg38_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg38_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh38_v0.6"

#' Somatic Hypermutation Locations GRCh38 Latest
#'
#' A data frame with somatic hypermutation locations in respect to GRCh38, the latest version.
#'
#' @format ## `somatic_hypermutation_locations_GRCh38_v_latest`
#' A data frame with 127 rows and 6 columns.
#' \describe{
#'   \item{chr_name}{Chromsome name.}
#'   \item{hg38_start}{Start coordinate for region in respect to hg19.}
#'   \item{hg38_end}{End coordinate for region in respect to hg19.}
#'   \item{gene}{Gene symbol in Hugo format.}
#'   \item{region}{Region name.}
#'   \item{regulatory_comment}{Annotates region with regulatory information.}
#' }
#' @keywords internal
"somatic_hypermutation_locations_GRCh38_v_latest"


#' GAMBL Metadata
#'
#' A data frame with metadata for a collection of GAMBL samples.
#'
#' @format ## `gambl_metadata`
#' A data frame with 4785 rows and 27 columns
#' \describe{
#'   \item{sample_id}{Sample identifier.}
#'   \item{patient_id}{Patient identifier.}
#'   \item{pathology}{Pathology.}
#'   \item{seq_type}{Sample sequencing type.}
#'   \item{genome_build}{Genome build the sample coordinates are in reference to.}
#'   \item{pairing_status}{Matched or unmatched.}
#'   \item{Tumor_Sample_Barcode}{Sample ID in another column, needed for certain functions.}
#'   \item{age_group}{Sample age group.}
#'   \item{compression}{The compression available for a particular sample.}
#'   \item{bam_available}{Boolean.}
#'   \item{pathology_rank}{Pathology rank.}
#'   \item{cohort}{Sample cohort}
#'   \item{COO_consensus}{COO consensus.}
#'   \item{DHITsig_consensus}{DHIT signature consensus.}
#'   \item{EBV_status_inf}{EBV status.}
#'   \item{ffpe_or_frozen}{FFPE or frozen.}
#'   \item{fl_grade}{FL grade.}
#'   \item{hiv_status}{Sample HIV status.}
#'   \item{lymphgen}{Lymphgen.}
#'   \item{lymphgen_cnv_noA53}{Lymphgen with CNV no A53.}
#'   \item{lymphgen_no_cnv}{Lymphgen no CNV.}
#'   \item{lymphgen_with_cnv}{Lymphgen with CNV}
#'   \item{lymphgen_wright}{Lymphgen Wright.}
#'   \item{molecular_BL}{Molecualr BL.}
#'   \item{normal_sample_id}{Normal sample ID}
#'   \item{sex}{Female or Male}
#'   \item{time_point}{Smaple timepoint.}
#' }
#' @keywords internal
"gambl_metadata"


#' Mirage Metrics
#'
#' A data frame with Quality Control metrics for a collection of GAMBL samples.
#' This resource is from the supplemental table of PMID: 36219881.
#'
#' @format ## `mirage_metrics`
#' A data frame with 2376 rows and 17 columns
#' \describe{
#'   \item{sample_id}{Sample identifier.}
#'   \item{AverageBaseQuality}{Average base quality.}
#'   \item{AverageInsertSize}{Average insert size.}
#'   \item{AverageReadLength}{Average read length.}
#'   \item{PairsOnDiffCHR}{Pairs on different chromosomes.}
#'   \item{TotalReads}{Total reads.}
#'   \item{TotallyUniquelyMapped}{Total uniquely mapped reads.}
#'   \item{TotalUnmappedreads}{Total unmapped reads.}
#'   \item{TotalDuplicatedreads}{Total duplicated reads.}
#'   \item{ProportionReadsDuplicated}{Proportion of duplicated reads.}
#'   \item{ProportionReadsMapped}{Proportion of mapepd reads.}
#'   \item{MeanCorrectedCoverage}{Mean corrected coverage.}
#'   \item{ProportionTargetsNoCoverage}{Proportion of targets with no coverage.}
#'   \item{ProportionCoverage10x}{Proportion of genome or target space with at least 10X coverage.}
#'   \item{study}{Study name}
#'   \item{coding_mutations}{Number of coding mutations.}
#'   \item{ProportionCoverage30x}{Proportion of genome or target space with at least 30X coverage.}
#' }
"mirage_metrics"


#' Hotspot Annotations
#'
#' A data frame with high-quality positions of ssm hotspots in selected genes.
#' This resource is based on GAMBL data and was used for hotspot annotation in
#' cFL/dFL classifier.
#'
#' @format ## `hotspots_annotations`
#' A data frame with 170 rows and 4 columns
#' \describe{
#'   \item{MAX_COORD}{Coordinate of hotspot region position identified by HotMaps.}
#'   \item{Start_Position}{Coordinate of gene position.}
#'   \item{Chromosome}{Name of gene chromosome.}
#'   \item{hot_spot}{Hot spot annotation.}
#' }
#' @keywords internal
"hotspots_annotations"


#' Protein Domains
#'
#' A data frame with high-quality positions of amino acid positions in their
#' corresponding domains.
#'
#' @format ## `protein_domains`
#' A data frame with 92849 rows and 12 columns
#' \describe{
#'   \item{HGNC}{HUGO Gene Nomenclature Committee.}
#'   \item{refseq.ID}{Reference sequence identifier.}
#'   \item{protein.ID}{Protein identifier.}
#'   \item{aa.length}{Amino acid length.}
#'   \item{Start}{Coordinate of amino acid position.}
#'   \item{End}{Coordinate of amino acid position.}
#'   \item{domain.source}{Domain source.}
#'   \item{Label}{Abbreviated domain type.}
#'   \item{domain.anno}{Domain type.}
#'   \item{pfam}{Protein family.}
#'   \item{Description}{Description of the domain.}
#'   \item{Description}{NA.}
#' }
#' @keywords internal
"protein_domains"

#' Cytobands coordinates (grch37)
#'
#' A data frame in bed format with coordinates of cytobands relative to grch37.
#'
#' @format ## `cytobands_grch37`
#' A data frame with 862 rows and 5 columns
#' \describe{
#'   \item{cb.chromosome}{Chromosome of the cytoband.}
#'   \item{cb.start}{Start position of the cytoband.}
#'   \item{cb.end}{End position of the cytoband.}
#'   \item{cb.name}{Cytoband name.}
#'   \item{label}{Cytoband label.}
#' }
#' @keywords internal
"cytobands_grch37"

#' Cytobands coordinates (hg38)
#'
#' A data frame in bed format with coordinates of cytobands relative to hg38.
#'
#' @format ## `cytobands_hg38`
#' A data frame with 862 rows and 5 columns
#' \describe{
#'   \item{cb.chromosome}{Chromosome of the cytoband.}
#'   \item{cb.start}{Start position of the cytoband.}
#'   \item{cb.end}{End position of the cytoband.}
#'   \item{cb.name}{Cytoband name.}
#'   \item{label}{Cytoband label.}
#' }
#' @keywords internal
"cytobands_hg38"

#' DLBCL90 genes
#'
#' A data frame with genes and their weights for DLBCL90.
#'
#' @format ## `dlbcl90_genes`
#' A data frame with 88 rows and 6 columns
#' \describe{
#'   \item{Gene.symbol}{Human-readable gene symbol matching the symbols used for the DLBCL90 code set design.}
#'   \item{coef}{Relative weight of the gene in DLBCL90 classification.}
#'   \item{Assay}{Name of the assay.}
#'   \item{ensembl_gene_id}{ENSEMBL gene id.}
#'   \item{gene_id}{ENSEMBL gene id with version.}
#'   \item{hgnc_symbol}{Human-readable gene symbol matching Gencode 33.}
#' }
#' @keywords internal
"dlbcl90_genes"

#' ENSG to symbol converter
#'
#' A data frame with ENSEMBL gene ids in different format and the corresponding
#'      hgnc symbols.
#'
#' @format ## `gencode_to_symbol`
#' A data frame with 60038 rows and 3 columns
#' \describe{
#'   \item{ensembl_gene_id}{ENSEMBL gene id in a format without .<number>.}
#'   \item{gene_id}{ENSEMBL gene id in a format containing .<number>.}
#'   \item{hgnc_symbol}{Human-readable gene symbol.}
#' }
#' @keywords internal
"gencode_to_symbol"

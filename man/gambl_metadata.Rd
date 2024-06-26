% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gambl_metadata}
\alias{gambl_metadata}
\title{GAMBL Metadata}
\format{
\subsection{\code{gambl_metadata}}{

A data frame with 4785 rows and 27 columns
\describe{
\item{sample_id}{Sample identifier.}
\item{patient_id}{Patient identifier.}
\item{pathology}{Pathology.}
\item{seq_type}{Sample sequencing type.}
\item{genome_build}{Genome build the sample coordinates are in reference to.}
\item{pairing_status}{Matched or unmatched.}
\item{Tumor_Sample_Barcode}{Sample ID in another column, needed for certain functions.}
\item{age_group}{Sample age group.}
\item{compression}{The compression available for a particular sample.}
\item{bam_available}{Boolean.}
\item{pathology_rank}{Pathology rank.}
\item{cohort}{Sample cohort}
\item{COO_consensus}{COO consensus.}
\item{DHITsig_consensus}{DHIT signature consensus.}
\item{EBV_status_inf}{EBV status.}
\item{ffpe_or_frozen}{FFPE or frozen.}
\item{fl_grade}{FL grade.}
\item{hiv_status}{Sample HIV status.}
\item{lymphgen}{Lymphgen.}
\item{lymphgen_cnv_noA53}{Lymphgen with CNV no A53.}
\item{lymphgen_no_cnv}{Lymphgen no CNV.}
\item{lymphgen_with_cnv}{Lymphgen with CNV}
\item{lymphgen_wright}{Lymphgen Wright.}
\item{molecular_BL}{Molecualr BL.}
\item{normal_sample_id}{Normal sample ID}
\item{sex}{Female or Male}
\item{time_point}{Smaple timepoint.}
}
}
}
\usage{
gambl_metadata
}
\description{
A data frame with metadata for a collection of GAMBL samples.
}
\keyword{datasets}

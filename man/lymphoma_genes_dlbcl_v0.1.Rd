% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{lymphoma_genes_dlbcl_v0.1}
\alias{lymphoma_genes_dlbcl_v0.1}
\title{Lymphoma Genes DLBCL v0.1}
\format{
\subsection{\code{lymphoma_genes_dlbcl_v0.1}}{

A data frame with 143 rows and 13 columns.
\describe{
\item{ensembl_gene_id}{Gene ID in ensembl format.}
\item{Gene}{Gene symbol in Hugo format.}
\item{Chappuy}{Boolean variable stating if the event is described in the study (Chappuy).}
\item{Reddy}{Boolean variable stating if the event is described in the study (Reddy).}
\item{LymphGen}{Boolean variable stating if the event is a described lymphgen or not.}
\item{curated}{Boolean variable annotating if the event is currated or not.}
\item{other_support}{Variable that annotates the event if there are other support available.}
\item{Lacy}{Boolean variable stating if the event is described in the study (Lacy).}
\item{aSHM}{Boolean varaible annotating if the event is considered an aSHM or not.}
\item{known_hotspots}{Boolean varaible annotating if the event is a known hotspot or not.}
\item{earliest_support}{Pubmeed ID to associated study.}
\item{common_alias}{Variable annotating other common aliases for the event, if such exists.}
\item{noncoding_driver_support}{Boolean variable annotating if the event has noncoding driver support or not.}
}
}
}
\usage{
lymphoma_genes_dlbcl_v0.1
}
\description{
Genes frequently associated with Diffuse large B cell lymphoma (DLBCL). This is version 0.1.
}
\keyword{internal}

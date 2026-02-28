#' Systemic Lupus Erythematosus (SLE) Gene Expression Data
#'
#' A subset of longitudinal gene expression data from a pediatric Systemic
#' Lupus Erythematosus (SLE) study. The full dataset contains 15,378 Illumina
#' HumanHT-12 V4.0 probes; this subset includes 500 probes plus 16 clinical
#' variables for a total of 519 columns. Loading this dataset creates an
#' object named \code{real_data}.
#'
#' @name SLE
#' @aliases real_data
#' @docType data
#' @format A data frame with 353 observations on 519 variables:
#' \describe{
#'   \item{id}{Subject ID (integer).}
#'   \item{y}{Response variable (continuous).}
#'   \item{intercept}{Intercept column (all ones).}
#'   \item{ILMN_*}{500 Illumina gene expression probes (numeric).}
#'   \item{AGE, WBC, NEUTROPHIL_COUNT, ESR}{Continuous clinical predictors.}
#'   \item{female, nonwhite}{Demographic indicators.}
#'   \item{ARTHRITIS, URINARY_CASTS, HEMATURIA, PROTEINURIA, PYURIA,
#'     NEW_RASH, MUCOSAL_ULCERS, LOW_COMPLEMENT, INCREASED_DNA_BINDING,
#'     LEUKOPENIA}{SLEDAI clinical components.}
#' }
#' @source Banchereau, R., Hong, S., Cantarel, B., et al. (2016).
#'   Personalized Immunomonitoring Uncovers Molecular Networks that Stratify
#'   Lupus Patients. \emph{Cell}, 165(3), 551--565.
#'   \doi{10.1016/j.cell.2016.05.057}.
#'   Gene Expression Omnibus accession
#'   \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391}{GSE65391}.
#' @usage data(SLE)
NULL

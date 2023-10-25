\name{SLE}
\alias{SLE}
\docType{data}
\title{
High-dimensional dataset for sparse linear mixed modeling.
}
\description{
This is the Systemic Lupus Erythematosus (SLE) data used in the data analysis section of the LMM-PROBE reference. The dataset has 309 observations, 103 subjects, and 15387 predictors. 
}
\usage{data(SLE)}
\format{
  A data frame with 309 observations and the following list elements:
  \describe{
    \item{\code{Y}}{First element of list, corresponding to the outcome to use in the lmmprobe function.}
    \item{\code{Z}}{Second element of list, corresponding to the high-dimesional matrix of sparse predictors for fixed effects.}
    \item{\code{V}}{Third element of list, corresponding to the low-dimesional matrix of non-sparse predictors for random effects. This matrix has either only one ID column, or one ID column with an additional column for a continuous variable for which a random slope is to be estimated.}
  }
}
\references{
Banchereau, R., Hong, S., Cantarel, B., Baldwin, N., Baisch, J., Edens, M., Cepika, A.-M., Acs, P., Turner, J., Anguiano, E., Vinod, P., Khan, S., Obermoser, G., Blankenship, D., Wakeland, E., Nassi, L., Gotte, A., Punaro, M., Liu, Y.-J., Banchereau, J., Rossello-Urgell, J., Wright, T., and Pascual, V. (2016), “Personalized Immunomonitoring Uncovers Molecular Networks that Stratify Lupus Patients,” Cell, 165, 551–565.

Zgodic, A. and McLain, A. C. (2023). Sparse high-dimensional linear mixed modeling with a partitioned empirical Bayes ECM algorithm. arXiv preprint arXiv:2310.12285.
}
\examples{
data(SLE)
Y = SLE$Y
Z = SLE$Z
V = SLE$V
}
\keyword{datasets}

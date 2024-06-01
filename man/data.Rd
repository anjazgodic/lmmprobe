\name{SLE}
\alias{SLE}
\docType{data}
\title{
High-dimensional dataset for sparse linear mixed modeling.
}
\description{
This is the Systemic Lupus Erythematosus (SLE) data used in the data analysis section of the LMM-PROBE reference. The dataset has 353 observations, 125 subjects, and 15397 predictors. 
}
\usage{data(SLE)}
\format{
  A data frame `real_data' with 353 observations.
}
\references{
Banchereau, R., Hong, S., Cantarel, B., Baldwin, N., Baisch, J., Edens, M., Cepika, A.-M., Acs, P., Turner, J., Anguiano, E., Vinod, P., Khan, S., Obermoser, G., Blankenship, D., Wakeland, E., Nassi, L., Gotte, A., Punaro, M., Liu, Y.-J., Banchereau, J., Rossello-Urgell, J., Wright, T., and Pascual, V. (2016), “Personalized Immunomonitoring Uncovers Molecular Networks that Stratify Lupus Patients,” Cell, 165, 551–565.

Zgodic A, Bay R, Zhang J, Olejua P, and McLain AC (2023). Sparse high-dimensional linear mixed modeling with a partitioned empirical Bayes ECM algorithm. arXiv preprint arXiv:2310.12285.
}
\examples{
data(SLE)
Y <- matrix(real_data[,"y"], ncol=1)
Z <- real_data[,4:ncol(real_data)]
V <- matrix(real_data[,"id"], ncol=1)
ID_data <- as.numeric(as.character(real_data$id))
}
\keyword{datasets}

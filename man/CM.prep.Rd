\name{CM.prep}
\alias{CM.prep}
\title{
Pre-process a gene expression dataset for CellMapper
}
\description{
Prepares a dataset for use with the \code{CM.search} function.
}
\usage{
CM.prep(Data, verbose = TRUE)
}
\arguments{
\item{Data}{
a numeric \code{matrix} or \code{ExpressionSet} object containing microarray expression data. The input expression data should be normalized and log-transformed.
}
\item{verbose}{
logical value indicating whether progress updates should be provided.
}
}
\details{
This function calculates the singular value decomposition of a gene expression matrix and prepares the data for use with the \code{CM.search} function. It can take some time for large expression matrices, and so it is recommended to save the output for future use. Pre-processed microarray data, ready for immediate use with \code{CM.search}, can be found in the \code{CellMapperData} package.

See the CellMapper vignette for examples about how to use \code{CM.prep} and associated functions to infer genes selectively expressed in specific cell types.
}
\value{
A list object ready to be provided to the \code{CM.search} function.
}
\references{
B.D. Nelms, L. Waldron, L.A. Barrera, A.W. Weflen, J.A. Goettel, G. Guo, R.K. Montgomery, M.R. Neutra, D.T. Breault, S.B. Snapper, S.H. Orkin, M.L. Bulyk, C. Huttenhower and W.I. Lencer. CellMapper: rapid and accurate inference of gene expression in difficult-to-isolate cell types. Genome Biol 2016, Accepted in Principle.
}
\author{
Brad Nelms
}
\seealso{
\code{\link{CM.search}}, \code{\link[Biobase]{ExpressionSet}}, \code{\link[CellMapperData]{CellMapperData}}
}
\examples{
# Create a mock expression dataset with random expression values
ngenes <- 1000
narrays <- 100
x <- matrix(rnorm(ngenes*narrays), ngenes, narrays)

# Prepare the dataset for use with CM.search
data <- CM.prep(x)

# Save the processed dataset for later
save(data, file = "Preprocessed_CellMapper_Data.RData")
}

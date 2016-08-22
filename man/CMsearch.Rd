\name{CMsearch}
\alias{CMsearch}
\alias{CellMapper}
\title{
Predict cell type-expressed genes using CellMapper
}
\description{
Predicts which genes are selectively expressed in the same cell type as a given 
cell type-specific marker gene (the 'query gene'), based on co-expression 
similarity.
}
\usage{
CMsearch(Data, query.genes = NULL, control.genes = NULL, QDW = TRUE, 
alpha = 0.5, verbose = TRUE, raw.pvals = FALSE)
}
\arguments{
\item{Data}{
a \code{CellMapperList} object containing expression data processed with 
\code{CMprep}, or a list of datasets processed with \code{CMprep}.
}
\item{query.genes}{
a list of genes that are specifically expressed in the cell type of interested, 
supplied as a character vector of gene names (matching the row names of the 
original expression matrix).
}
\item{control.genes}{
a list of genes expressed specifically in non-target cell types (optional), 
supplied as a character vector of gene names (matching the row names of the 
original expression matrix). This option generally has little effect on the 
results and its use is NOT recommended.
}
\item{QDW}{
logical value indicating whether 'query driven weighting' should be applied in 
the CellMapper SVD filter. The default value of TRUE can be used in most cases.
}
\item{alpha}{
alpha parameter controlling the strength of the CellMapper SVD filter. The 
default value of 0.5 can be used in most cases.
}
\item{verbose}{
logical value indicating whether progress updates should be provided.
}
\item{raw.pvals}{
logical value indicating whether unadjusted p-values, which have not been 
corrected for multiple hypothesis testing, should be returned.
}
}
\details{
This function predicts which genes are selectively expressed in the same cell 
type as a given cell type-specific marker gene (the 'query gene'), based on 
co-expression similarity. The only required inputs are a gene expression matrix 
that has been pre-processed with the \code{CMprep} function (or a list of 
pre-processed expression matrices), and a 'query gene' known to be specifically 
expressed in the cell type of interest. Pre-processed microarray data, ready for
 immediate use with \code{CMsearch}, can be found in the \code{CellMapperData} 
package.

See the CellMapper vignette for examples about how to use \code{CMsearch} and 
associated functions to infer genes selectively expressed in specific cell 
types.
}
\value{
A dataframe containing the gene identifiers in the first column, the false 
discovery rate (FDR) in the second, and the unadjusted p-value in the third. 
FDR is calculated using the method of Benjamini and Hochberg.
}
\references{
B.D. Nelms, L. Waldron, L.A. Barrera, A.W. Weflen, J.A. Goettel, G. Guo, R.K. 
Montgomery, M.R. Neutra, D.T. Breault, S.B. Snapper, S.H. Orkin, M.L. Bulyk, C. 
Huttenhower and W.I. Lencer. CellMapper: rapid and accurate inference of gene 
expression in difficult-to-isolate cell types. Genome Biol 2016, Accepted in 
Principle.
}
\author{
Brad Nelms
}
\seealso{
\code{\link{CMprep}}, \code{\link[CellMapperData]{CellMapperData}}
}
\examples{
# Create a mock expression dataset with random expression values
ngenes <- 1000
narrays <- 100
x <- matrix(rnorm(ngenes*narrays), ngenes, narrays)
rownames(x) <- 1:ngenes

# Prepare the dataset for use with CMsearch
data <- CMprep(x)

# Predict which genes are co-expressed in the same cell type as 'gene' 59 in the
# mock expression dataset
results <- CMsearch(data, query.genes = '59')

head(results)
}

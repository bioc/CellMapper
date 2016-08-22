\name{CellMapperList-class}
\docType{class}

\alias{CellMapperList}
\alias{CellMapperList-class}
\alias{show,CellMapperList-method}

\title{
Class to contain expression data pre-processed for the 'CellMapper' package. The
 'CellMapperList' class extends the 'SimpleList' class.
}
\description{
Container for expression data that has been pre-processed using the 'CMprep' 
function. A 'CellMapperList' object can be provided directly to the 'CMsearch' 
function to predict genes expressed selectively in specific cell types.
}
\details{
'CMprep' transforms an expression matrix using singular value decomposition 
(SVD), resulting in a matrix, 'B', with the left-singular vectors of original 
data matrix and a vector, 'd', with the singular values. Singular vectors that 
account for less variance than an individual sample in the original dataset have
 been trimed (Kaiser's criterion), thereby removing singular vectors that mainly
 account for noise. The advantage of this transformation is that it reduces 
dataset size, and avoids the need to perform a time-consuming SVD transformation
 before running 'CMsearch'.

A 'CellMapperList' object contains the transformed left singular value matrix, 
'B', and singular values, 'd', as well as meta-data associated with the original
 expression matrix.
}
\section{Constructor:}{
'CellMapperList' instances are usually created through the 'CMprep' function. 
See ?CMprep

To create a 'CellMapperList' object directly, the following constructor is 
provided:

CellMapperList(B, d, meta = list())

where B is a numeric matrix containing left-singular vectors, d is a numeric 
vector containing singular values, and meta is a list object containing 
metadata.
}
\section{Accessors}{
Same as for SimpleList objects. See ?SimpleList

The sample metadata can be extracted using metadata().
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
\code{\link{CMsearch}}, \code{\link{CMprep}}
}
\examples{
# Create a mock expression dataset with random expression values
ngenes <- 1000
narrays <- 100
x <- matrix(rnorm(ngenes*narrays), ngenes, narrays)
rownames(x) <- 1:ngenes

# Prepare a CellMapperList object using the CMprep function
data <- CMprep(x, DataSource = "Mock Expression Matrix")
show(data)
metadata(data)
}
\keyword{classes}

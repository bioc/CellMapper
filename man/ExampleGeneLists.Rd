\name{ExampleGeneLists}
\alias{ExampleGeneLists}
\alias{ExampleQueryGenes}
\alias{PositiveControlGenes}
\alias{NegativeControlGenes}
\docType{data}
\title{
Example Gene Lists
}
\description{
Query genes, positive control genes, and negative control genes for example cell
 types.
}
\usage{
data(ExampleQueryGenes)
data(PositiveControlGenes)
data(NegativeControlGenes)
}
\details{
'ExampleQueryGenes' is a data frame listing the cell types and query genes 
analyzed in Nelms, et al. (2016). These genes can be used as query or control 
genes for custom CellMapper searches, and are also provided for reproducibility 
of the manuscript results. All predictions from the manuscript can be generated 
using the \code{ReplicateManuscript} function. The 'Dataset' column indicates 
which dataset(s) were used for each cell type, available from Bioconductor's 
\code{ExperimentHub}. All datasets are available from the \code{CellMapperData} 
package.

'PositiveControlGenes' is a list of positive control genes for each cell type 
used in Nelms, et al. (2016). These positive control genes were NOT provided to 
the algorithm, but were used to assess the accuracy of the results. This list 
contains Entrez gene IDs.

'NegativeControlGenes' is a list of negative control genes for each cell type 
used in Nelms, et al. (2016). These negative control genes were NOT provided to 
the algorithm, but were used to assess the accuracy of the results. This list 
contains Entrez gene IDs.
}
\value{
'ExampleQueryGenes' is a data frame, 'PositiveControlGenes' and 
'NegativeControlGenes' are both lists of character vectors
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
\code{\link{CMsearch}}, \code{\link{ReplicateManuscript}}, 
\code{\link[CellMapperData]{CellMapperData}}, 
\code{\link[ExperimentHub]{ExperimentHub}}
}
\examples{
data(ExampleQueryGenes)
head(ExampleQueryGenes)
}

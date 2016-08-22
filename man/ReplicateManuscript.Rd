\name{ReplicateManuscript}
\alias{ReplicateManuscript}
\title{
Replicate Results from CellMapper Manuscript
}
\description{
Recreates the CellMapper predictions from Nelms, et al. (2016) -- Additional 
dataset 1.
}
\usage{
ReplicateManuscript()
}
\details{
This function will recreate the CellMapper predictions from Nelms, et al. 
(2016), Additional dataset 1. It requires the \code{ExperimentHub} and 
\code{CellMapperData} packages to be installed. Provided that the 
\code{CellMapperData} resources have previously been downloaded from 
\code{ExperimentHub}, the 'ReplicateManuscript' function should take ~5 minutes 
to run on a personal laptop.
}
\value{
A 'list' object containing CellMapper predictions for each cell type. Genes 
names are returned as Entrez identifiers.
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
\code{\link{CMsearch}}, \code{\link{ExampleQueryGenes}}, 
\code{\link[CellMapperData]{CellMapperData}}, 
\code{\link[ExperimentHub]{ExperimentHub}}
}
\examples{
Results <- ReplicateManuscript()
head(Results[['GABAergic Neurons']])
}

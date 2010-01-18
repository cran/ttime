\name{phylo}
\Rdversion{1.1}
\alias{phylo}
\title{
 Phylogenetic proximity from neurodevelopmental event timings.
}

\description{
Hierarchical clustering with euclidean metric and complete linkage establishing possible phylogenetic proximity between the species from the known and predicted events (Nagarajan and Clancy, 2008). Uncertainty in the clustering results is investigated using pvclust. 
}

\usage{
  phylo(out)
}

\arguments{
  \item{out}{
	Predicted neurdevelopmental events.
   }
}

\value{
Dendrogram establishing possible phylogenetic proximity between the species based on hierachical clustering of the known and predicted event timings. Complete linkage and euclidean metric are used are arguements in constructing the dendrogram.
}


\references{Nagarajan, R and Clancy, B (2008) \emph{Phylogenetic Proximity Revealed by Neurodevelopmental Event Timings.} Neuroinformatics, 6(2), 71-79. \cr
Suzuki, R Shimodaira, H (2006) \emph{Pvclust: an R package for assessing the uncertainty in hierarchical clustering.} Bioinformatics 22, 1540-1542.

}

\author{
Radhakrishnan Nagarajan
}

\seealso{\code{\link{event_data}}, \code{\link{translate}}}

\examples{
library(ttime);
data(event_data); # event timing data
npsp <- 8;  # number of non-primate species
out <- translate(event_data, npsp); # predicted events
phylo(out); # hierarchical clustering
}

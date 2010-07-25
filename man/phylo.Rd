\name{phylo}
\Rdversion{1.2.1}
\alias{phylo}
\title{
 Phylogenetic proximity between species inferred from their neurodevelopmental event timing
}

\description{

Dendrogram establishing possible phylogenetic proximity between the species based on hierachical clustering of the known and predicted event timings  (Nagarajan and Clancy 2008) returned by the function \emph{translate}. Complete linkage and euclidean metric are used for generating the dendrogram. Uncertainty in the clustering results is reflected by the \emph{au} (approximately unbiased) and \emph{bp} (bootstrap) probabilities (Shimodaira and Suzuki 2006).

}

\usage{
  phylo(data)
}

\arguments{
  \item{data}{
	A data frame containing the known and predicted neurodevelopmental event timing across species returned by the function \emph{translate}. 
   }
}

\value{
The function generates a dendrogram and does return any values (NULL). 
}


\references{
Nagarajan, R., & Clancy, B. (2008). \emph{Phylogenetic proximity revealed by neurodevelopmental event timings.} Neuroinformatics, 6, 71-79.\cr

Suzuki, R., Shimodaira, H. (2006). \emph{Pvclust: an R package for assessing the uncertainty in hierarchical clustering.} Bioinformatics, 22,1540-1542.\cr

}

\author{
Radhakrishnan Nagarajan
}

\seealso{\code{\link{event_data}}, \code{\link{translate}}}

\examples{
library(ttime);
data(event_data); # event timing data
npsp <- 8;  	  # number of non-primate species
data <- translate(event_data, npsp); # predicted events
phylo(data); 	  # hierarchical clustering
}

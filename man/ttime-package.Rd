\name{ttime-package}
\Rdversion{1.2.4}
\alias{ttime-package}
\alias{ttime}
\docType{package}
\title{
Translate Neurodevelopmental Event Timing Across Species
}
\description{
Neurdevelopmental event timing sequence is conserved across species (Finlay and Darlington 1995). However, only a handful of event timings have been determined experimentally. The ttime package implements the necessary functions to translate neurdoevelopmental event timing across species (Finlay and Darlington 1995; Clancy et al., 2000) and investigate possible phylogenetic relationships (Nagarajan and Clancy 2008). Unknown events across species are predicted using the model \emph{score + event score = ln(post-conceptional day - k)} proposed in (Finlay and Darlington 1995; Clancy et al., 2000).

}

\details{
\tabular{ll}{
Package: \tab ttime\cr
Type: \tab Package\cr
Version: \tab 1.2.4\cr
Date: \tab 2010-08-03\cr
License: \tab GPL >=2.0 \cr
LazyLoad: \tab yes
 }
}

\author{
Radhakrishnan Nagarajan\cr
Maintainer: <rnagarajan@uams.edu>
}
\references{
Finlay, B.L., & Darlington, R.B. (1995). \emph{Linked regularities in the development and evolution of mammalian brains.} Science, 268,1578-1584.\cr

Clancy, B., Darlington, R.B., & Finlay, B.L. (2000).  \emph{The course of human events:  predicting the timing of primate neural development.} Developmental Science, 3, 57-66.\cr

Nagarajan, R., & Clancy, B. (2008). \emph{Phylogenetic proximity revealed by neurodevelopmental event timings.} Neuroinformatics, 6, 71-79.\cr

Suzuki, R., Shimodaira, H. (2006). \emph{Pvclust: an R package for assessing the uncertainty in hierarchical clustering.} Bioinformatics, 22,1540-1542.\cr
}

\seealso{\code{\link{event_data}}, \code{\link{translate}}, \code{\link{phylo}}}



\examples{
library(ttime);
#Event timing data consisting of known and unknown events 
data(event_data); 
#Number of non-primate species in event_data
npsp <- 8;   
#Predict the unknown event timings and their 95% confidence interval
pred_vals <- translate(event_data, npsp);  
#Dendrogram by hierachical clustering of the event timings
phylo(pred_vals); 
}

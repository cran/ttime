\name{ttime-package}
\Rdversion{1.1}
\alias{ttime-package}
\alias{ttime}
\docType{package}
\title{
Translate Neurodevelopmental Event Timing Across Species
}
\description{
Translate neurdoevelopmental event timing across species and predict unknown neurodevelopmental event timing based on the model  species score + event score = ln(post-conceptional day - k) proposed in (Finlay and Darlington, 1995; Clancy et al., 2000). Investigate phylogenetic proximity between the species from the known and predicted event timings (Nagarajan and Clancy, 2008).
}
\details{
\tabular{ll}{
Package: \tab ttime\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-01-15\cr
License: \tab GPL >=2.0 \cr
LazyLoad: \tab yes\cr
}

}
\author{
Radhakrishnan Nagarajan\cr
Maintainer: <rnagarajan@uams.edu>
}
\references{Finlay BL, Darlington RB. (1995) \emph{Linked regularities in the development and evolution of mammalian brains.} Science, 268 (5217): 1578-84.\cr
Clancy B, Darlington RB, Finlay BL. (2000) \emph{The course of human events:  predicting the timing of primate neural development.} Developmental Science 2000;3:57-66.\cr
Nagarajan, R and Clancy, B (2008) \emph{Phylogenetic Proximity Revealed by Neurodevelopmental Event Timings.} Neuroinformatics, 6(2), 71-79.
}

\seealso{\code{\link{event_data}}, \code{\link{translate}}, \code{\link{phylo}}}



\examples{
library(ttime);
#Event timing data consisting of known and unknown events
data(event_data); 
#Number of non-primate species in event_data
npsp <- 8;   
#Predict the unknown event timings with 95\% confidence interval
out <- translate(event_data, npsp);  
#Hierachical clustering of the known and predicted event timings
phylo(out); 
}

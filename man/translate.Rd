\name{translate}
\Rdversion{1.1}
\alias{translate}
\title{
Translates neurodevelopmental event timing across species and predicts unknown event timings.
}
\description{
Neurdevelopmental event timing sequence is conserved across species (Finlay and Darlington, 1995). However, only a handful of event timings have been determined experimentally. The function translate predicts the unknown events using the model species score + event score = ln(post-conceptional day - k) proposed in (Finlay and Darlington, 1995; Clancy et al., 2000).
}
\usage{
out <- translate(yo, npsp)
}
\arguments{
  \item{yo}{
	A data frame containing the neurodevelopmental event timing across species. 
}
  \item{npsp}{
       Number of non-primate species in yo.
}
}
\details{
A sufficient number of known events (non-zero) are necessary to obtain reliable predictions.
The structure of yo should be the same as that of event_data.
}

\references{Finlay BL, Darlington RB. (1995) \emph{Linked regularities in the development and evolution of mammalian brains.} Science, 268 (5217): 1578-84.\cr
Clancy B, Darlington RB, Finlay BL. (2000) \emph{The course of human events:  predicting the timing of primate neural development.} Developmental Science 2000;3:57-66.}



\author{
Radhakrishnan Nagarajan
}

\seealso{\code{\link{event_data}}}

\examples{
library(ttime);
data(event_data); #event timing data
npsp <- 8;  # number of non-primate species in event_data
out <- translate(event_data, npsp); #predicted events with 95\% confidence interval
}

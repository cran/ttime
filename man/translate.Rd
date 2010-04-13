\name{translate}
\Rdversion{1.1}
\alias{translate}
\title{
Translate neurodevelopmental event timing across species
}
\description{
The function translate predicts the unknown events using the model \emph{species score + event score = ln(post-conceptional day - k)} proposed in (Finlay and Darlington, 1995; Clancy et al., 2000). The parameter \emph{k} is estimated from the data by maximizing the correlation between the known and predicted events (Clancy et al., 2000).
}
\usage{
pred_vals <- translate(inp_vals, nonprim)
}
\arguments{
  \item{inp_vals}{A data frame containing the known and unknown neurodevelopmental event timing across species. 
}
  \item{nonprim}{Number of non-primate species in \emph{inp_vals}.
}
}

\value{The known and predicted events along with their 95\% confidence interval are returned in \emph{pred_vals}.

}
\details{
A sufficient number of known events are necessary to obtain reliable predictions. The structure of \emph{input_vals} should be the same as that of event_data.
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
nonprim <- 8;  # number of non-primate species in event_data
out <- translate(event_data, nonprim); #predicted events with 95\% confidence interval
}

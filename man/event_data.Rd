\name{event_data}
\Rdversion{1.1}
\alias{event_data}
\docType{data}
\title{
Neurodevelopmental event timing data consisting of 106 known and unknown events across 10 mammals (8 non-primates and 2 primates).
}
\description{
Neurodevelopmental event timing data consisting of 106 known and unknown events across 10 mammals. Unknown events are represented by zeros.
}
\usage{data(event_data)}
\format{
  Variables in event_data follow the order below.
  \describe{
    \item{\code{Event}}{Name of the neurodevelopmental event}
    \item{\code{Hamster}}{Hamster neurodevelopmental event timing}
    \item{\code{Mouse}}{Mouse neurodevelopmental event timing}
    \item{\code{Rat}}{Rat neurodevelopmental event timing}
    \item{\code{Rabbit}}{Rabbit neurodevelopmental event timing}
    \item{\code{Spiny.Mouse}}{Spiny-Mouse neurodevelopmental event timing}
    \item{\code{Guinea.Pig}}{Guinea Pig neurodevelopmental event timing}
    \item{\code{Ferret}}{Ferret neurodevelopmental event timing}
    \item{\code{Cat}}{Cat neurodevelopmental event timing}
    \item{\code{Macaque}}{Macaque neurodevelopmental event timing}
    \item{\code{Human}}{Human neurodevelopmental event timing}
    \item{\code{Cortical}}{Cortical neurodevelopmental event (1 yes, 0 no) }
    \item{\code{Limbic}}{Limbic neurodevelopmental event (1 yes, 0 no) }
    \item{\code{Reference}}{Bibiographic Reference for the corresponding event}
 }
}
\examples{
#Neurodevelopmental event timing data consisting of 106 known and unknown events across 10 species.
#The columns are arranged in the order described above
library(ttime);
data(event_data);
}
\seealso{http://www.translatingtime.net/}

\references{Finlay BL, Darlington RB. (1995) \emph{Linked regularities in the development and evolution of mammalian brains.} Science, 268 (5217): 1578-84.\cr
Clancy B, Darlington RB, Finlay BL. (2000) \emph{The course of human events:  predicting the timing of primate neural development.} Developmental Science 2000;3:57-66.}



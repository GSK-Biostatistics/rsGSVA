% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_data.R
\name{simulate_data}
\alias{simulate_data}
\title{Simulate data with gene sets}
\usage{
simulate_data(
  pPos = 1,
  n = 100,
  pfh = 10,
  pnh = 10,
  mnh = 5,
  pfi = 10,
  pni = 100,
  SNRh = 0.5,
  SNRo = 0.5,
  rFHub = 0.3,
  rNFHub = 0.3,
  sdNoise = 1,
  vrb = F
)
}
\arguments{
\item{pPos}{Probability of positive coefficient for indicator variables}

\item{n}{Total number of samples.}

\item{pfh}{The number of features in the functional hub}

\item{pnh}{The number of features in each non-functional hub.}

\item{mnh}{The number of non-functional hubs}

\item{pfi}{The number of features that are functional independents.}

\item{pni}{The number of features that are nonfunctional independents.
Note: The total number of features will be denoted by p where
p=pfh+mnh*pnh+pfi+pni.}

\item{SNRh}{Hub signal to noise ratio. Number between 0 and 1  (1 indicates high signal)}

\item{SNRo}{Signal to noise ratio for a functional independent. Number between 0 and 1 (1 indicates high signal)}

\item{rFHub}{Correlation parameter between features within the functional hub}

\item{rNFHub}{Correlation between features corresponding to the nonfunctional
hub and each of their respective hubs (see xNFHub below)
note: for the time being, betaNFHub is the same for all
non functional hubsalpha = correlation parameter between xhub values (see below) and Y.}

\item{sdNoise}{Standard deviation for the non-functional independents}

\item{vrb}{Verbose flag to enable printing of intermediate steps}
}
\value{
{ object of class "simHub" containing
x = p x n data matrix (see below).
y = n x 1 vector.
featureID=px1 character vector with values ('fh','nh','fi','ni').

NOTE: the rows of x are ordered as follows:

functional hub features      1:pfh

non-functional hub features  (pfh+1):((pfh+1)+mnh*pnh)

functional independents      ((pfh+1)+mnh\emph{pnh+1):(((pfh+1)+mnh}pnh+1)+pfi)

non-functional independents  (((pfh+1)+mnh*pnh+1)+pfi+1):p

standard deviation of background noise in y is 1.0}
}
\description{
\code{simulate_data} simulates data with gene sets for different relationships. A detailed description is given in Miecznikowski et.al.
}
\details{
Simulates the following gene sets:
\itemize{
\item Functional hub (genes prefixed with "fh"): Genes that are correlated with each other and the outcome
\item Non-functional hubs (genes prefixed with "nh"): Genes that are correlated with each other but not the outcome
\item Functional independents (genes prefixed with "fi"): Genes that are correlated with the outcome but not each other
\item Non-functional independents (genes prefixed with "ni"): Genes that are not correlated with each other or the outcome. Independent noise.
}
}
\references{
Miecznikowski JC, Gaile DP, Chen X, Tritchler DL.
Identification of consistent functional genetic modules.
\emph{Stat Appl Genet Mol Biol.} 2016 Mar;15(1):1-18. doi:
10.1515/sagmb-2015-0026.
}
\author{
Lorin Towle-Miller
}

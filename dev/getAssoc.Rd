\name{getAssoc}
\alias{getAssoc}
\title{\code{getAssoc} function}
\usage{
  getAssoc(res,ATC=NULL,meddra=NULL)
}
\arguments{
  \item{res}{must be the output of one the signal detection
  method function GPS.}

  \item{ATC}{The label of the drug. By default, ATC=FALSE.}

  \item{meddra}{The label of the adverse event. By default,
  meddra=FALSE.}
}
\value{
  ATC Recalls the label of the drug.

  meddra Recalls the label of the event.

  EXIST_ATC Indicates if the label of the drug exists in
  the database.

  EXIST_meddra Indicates if the label of the adverse event
  exists in the database.

  EXIST_COUPLE Indicates if the couple is present in the
  database.

  LIST It is a dataframe that contains the labels, the
  counts, the expected counts, the value of the statistic
  of interest, the rank and the estimated FDR for each
  couple.
}
\description{
  This function makes possible to extract some information
  from the output of the ADRminer functions for a given
  couple adverse event-drug, for a drug or for an adverse
  event.
}
\author{
  Youness Ergaibi & Ismaïl Ahmed & Antoine Poncet
}
\keyword{getAssoc}


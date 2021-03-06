% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpectroX.R
\name{getPeptides}
\alias{getPeptides}
\title{Get Mz Shift of complementary peptide/fragment ion}
\usage{
getPeptides(proteinSeq, proteaseRegExp = "[KR](?!P)", nbMiscleavages = 0)
}
\arguments{
\item{proteinSeq}{protein sequence}

\item{proteaseRegExp}{protease Regular Expression}

\item{nbMiscleavages}{default 0}

\item{aaSeq}{character string}

\item{charge}{default 1}

\item{isHeavy}{c(T,F)}

\item{annotSpectrum}{Annotated Spectrum object}
}
\value{
numeric

annotSpectrum Annotated Spectrum object

vector of peptides
}
\description{
Get Mz Shift of complementary peptide/fragment ion
}
\details{
Calculate mass shift of complimentary spectrum.
 \itemize{
 \item PETIDEK (light) -> 8.014199
 \item PETIDEK (heavy) -> -8.014199
 \item PETIDER (light) -> 10.008269
 \item PETIDER (heavy) -> -10.008269
 }

No details

No details
}
\note{
No note

No note

No note
}
\examples{
print("No examples")
Created complimentary Heavy/Light annotated spectrum with updates precursor and fragment m/z values.
print("No examples")
Digest protein
print("No examples")
}
\references{
NA

NA
}

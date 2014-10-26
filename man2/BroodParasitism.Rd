\name{BroodParasitism}
\alias{BroodParasitism}
\docType{data}
\title{
	Conspecific Brood Parasitism Dataset
}
\description{
	Dataset containing information on intraspecific brood parasitsim in a population of domesticated zebra finches.
}
\usage{
	data(BroodParasitism) 
} 
\format{
  A data frame with 98 observations on the following 9 variables.
  \describe{
    \item{FemaleID}{Female bird identities (49 individuals in total).}
    \item{Season}{One out of two sesaons (2005 and 2006). Note that some individuals were present in only one season.}
    \item{PairedYN}{Whether or not a particular female was socially monogomously paired in a given season.}
	\item{PairedBothSeasons}{Whether or not a female was paird in both seasons.}
	\item{OwnClutchesBothSeasons}{Whether or not a female was incubating own clutches in both season.}
    \item{OwnClutches}{The number of own clutches produced and incubated with a social partner.}
    \item{HostClutches}{The number of own clutches that contained at least one parasitic egg (laid by other females).}
    \item{HostYN}{Binary score for being a host in a given season. Is =1 if \code{ParasitisedClutches>0} and =0 if \code{ParasitisedClutches==0}.}
    \item{cbpEggs}{Number of parasitic eggs dumped to clutches of other females.}
    \item{cbpYN}{Binary score of whether or not a female layed parasitic eggs in a given season. Is =1 if \code{dumpedEggs>0} and =0 if \code{dumpedEggs=0}.}
    \item{nEggs}{Total number of own eggs laid in a given season.}
  }
}
\details{
	Females were allowed to breed for three month in 2005 and for three month in 2006. 
	They lived in aviaries with 6-9 females and 6-9 males per aviary. 
	Individuals were swapped among aviaries between the 2005 and the 2006 breeding season, but some females were used in one season only. 
	Hatchling and embryos were genotyped and parentage was assigned by exclusion. 
	About 21% of all clutches that were initiated and incubated by a social pair contained at least one dumped egg (that was laid by another female).
	Many females applied a mixed strategy combining caring for own clutches with brood parasitism.
	Some unpaired (or part-time paired) females also applied brood parasitism.
}
\references{
	Schielzeth, H. and Bolund E. (2010) \emph{Patterns of conspecific brood parasitism in zebra finches}. Animal Behaviour 79: 1329-1337.
} 
\keyword{datasets}
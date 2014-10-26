\name{print.rpt}
\alias{print.rpt}
\title{
	Print a rpt object
}
\description{
	Displays the results a rpt object (i.e. the result of a rpt function call) an a nice form.
}
\usage{
	\method{print}{rpt}(x, \dots)
}
\arguments{
  \item{x}{An rpt object returned from one of the rpt functions}
  \item{\dots}{Additional arguments; none are used in this method.}
}
\references{
	Nakagawa, S. and Schielzeth, H. (2011) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956.
}
\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\seealso{
	\link{rpt}, \link{rpt.corr}, \link{rpt.aov}, \link{rpt.remlLMM}, \link{rpt.mcmcLMM}, \link{rpt.binomGLMM.add}, \link{rpt.binomGLMM.multi}, \link{rpt.poisGLMM.add}, \link{rpt.poisGLMM.multi} 
}
\examples{
# repeatability estimation for weight (body mass)
data(BodySize)
attach(BodySize)
print.rpt(rpt.Weight <- rpt.mcmcLMM(Weight, BirdID))
print(rpt.Weight)  # alternative call to print.rpt() through pring()
detach(BodySize)
}
\keyword{models}
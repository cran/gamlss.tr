\name{gamlss.tr-package}
\alias{gamlss.tr-package}
\alias{gamlss.tr}
\docType{package}
\title{The GAMLSS add on package for truncated distributions}
\description{
The main purpose of this package is to allow the user of the GAMLSS models to fit truncated distributions.
}
\details{
\tabular{ll}{
Package: \tab gamlss.tr\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2005-11-21\cr
License: \tab GPL (version 2 or later) \cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Mikis Stasinopoulos <\email{d.stasinopoulos@londonmet.ac.uk}> and Bob Rigby <\email{r.rigby@londonmet.ac.uk}> 

Maintainer: Mikis Stasinopoulos <d.stasinopoulos@londonmet.ac.uk>

}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.londonmet.ac.uk/gamlss/}).
}
\keyword{ package }
\keyword{ distribution }
\keyword{ regression }

\seealso{\code{\link[gamlss]{gamlss}},\code{\link[gamlss]{gamlss.family}}  
}
\examples{
gen.trun(par=c(0,100),family="TF", name="0to100", type="both")
plot(function(x) dTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100)
plot(function(x) pTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100)
plot(function(x) qTF0to100(x, mu=80 ,sigma=20, nu=5), 0.01, .999)
hist(rTF0to100(1000, mu=80 ,sigma=20, nu=5))
}

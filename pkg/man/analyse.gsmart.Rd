\name{analyse.gsmart}
\alias{analyse.gsmart}
\title{
Main Function of the gsmart package.
}
\description{
Allows for mutliple gene set testing despite inter-gene and inter-gene set correlations as proposed in Tian et al.

}
\usage{
analyse.gsmart(X, response, gs.list, nrot = 1000, tests = c("W", "G"), one.sided = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Expression matrix X
%%     ~~Describe \code{X} here~~
}
  \item{response}{
Response Vector
%%     ~~Describe \code{response} here~~
}
  \item{gs.list}{
Data frame or matrix of gene sets. Genes in first, gene sets in second column.
%%     ~~Describe \code{gs.list} here~~
}
  \item{nrot}{
Number of rotations applied.
%%     ~~Describe \code{nrot} here~~
}
  \item{tests}{
Which tests shall be applied? Currently, Wilcoxon ("W") or Goeman ("G") can be chosen. Put "all" if you want all tests to be performed.
%%     ~~Describe \code{tests} here~~
}
  \item{one.sided}{
logical, optional. Shall one-sided p-values be returned as well? Only available for Wilcoxon.
%%     ~~Describe \code{one.sided} here~~
}
  \item{verbose}{
Extra output.
%%     ~~Describe \code{verbose} here~~
}
}

\value{
result list
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Mathias Fuchs and Stephan Artmann, in preparation
Lu Tian, Steven A. Greenberg, Sek Won Kong, Josiah Altschuler, Isaac S. Kohane, and Peter J. Park Discovering statistically significant pathways in expression profiling studies

%% ~put references to the literature/web site here ~
}
\author{
Stephan Artmann
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function (X, response, gs.list, nrot = 1000, tests = c("W", "G"), 
#    one.sided = FALSE, verbose = FALSE) 
#{
#    R = test.gs(X, response, gs.list, nrot = nrot, tests = tests, 
#        one.sided = one.sided, verbose = verbose)
#    Res = list()
#    for (i in 1:length(names(R$S.star))) {
#        Res[[i]] = do.MHC.R(R$S.star[[i]])
#    }
#    names(Res) = names(R$S.star)%    Res
#  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

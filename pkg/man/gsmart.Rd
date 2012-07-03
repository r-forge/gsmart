\name{gsmart}
\alias{gsmart}
\docType{package}
\title{
Package Description: Gene Set MicroArray Testing.
}
\description{
Test mRNA expression values against an alternative when mRNAs are grouped into gene sets. 
}
\details{
\tabular{ll}{
Package: \tab gsmart\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2012-03-05\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
For a detailed help check vignette("gsmart")

You can start the test with the "analyse.gsmart" which has its own help file.
}
\author{
Stephan Artmann <stephan-artmann@gmx.net>,
Mathias Fuchs

Maintainer: Stephan Artmann <stephan-artmann@gmx.net>
}
\references{
Mathias Fuchs and Stephan Artmann, Gene Set MicroArray Testing (in preparation).
}
\keyword{ package, gene sets
mRNA testing,
expression }
\seealso{
Function "analyse.gsmart"
}
\examples{
#######################################
### Generate random expression data ###
#######################################
# Generate random mRNA expression data with 20 mRNAs
# and 10 replicates
Y = rnorm(200);
dim(Y) = c(20,10);
rownames(Y) = 1:20;

# Lets assume that the corresponding mRNA experiments had 5 replicates in each group
groups = (c(1,1,1,1,1,2,2,2,2,2));

####################
### Perform Test ###
####################
library(gsmart)
#Let miRNA 1 attack mRNAs 1 to 9 and miRNA 2 attack mRNAs 10 to 17.
# mRNAs 18 to 20 are not attacked. miRNA 3 has no gene set.
miR = c(rep(1,9),c(rep(2,8)));
mRNAs = 1:17;
A = data.frame(mRNAs,miR); # Note that the gene sets MUST be in the second column!
A

###############
### Analyse ###
###############
set.seed(1)
P = analyse.gsmart(Y,groups,A)
P

}






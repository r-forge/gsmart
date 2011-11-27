
#' Main Function
#' @param X Expression matrix X.
#' @param response Response vector.
#' @param gs.list Data frame or matrix of gene sets. Genes in first, gene sets in second column.
#' @param nrot Number of rotations to be performed
#' @param tests Which tests shall be applied? Currently, Wilcoxon ("W") or Goeman ("G") can be chosen. Put "all" if you want all tests to be performed.
#' @param one.sided logical, optional. Shall one-sided p-values be returned as well? Only available for Wilcoxon.
#' @return Tables of results.
#' @author Stephan Artmann
#' @examples
##MAINEXAMPLE
analyse.gs = function(X,response,gs.list,nrot=1000,tests=c("W","G"), one.sided = FALSE, verbose = FALSE) {
 R = test.gene.sets(X,response,gs.list,nrot=nrot,tests=tests, one.sided = one.sided, verbose = verbose);
 Res = list();
 for (i in 1:length(names(R$S.star))) {
  Res [[i]] = do.MHC.R(R$S.star[[i]]);
 }
 names(Res) = names(R$S.star);
 Res;
}

#' Generation of Result Table.
#' @param S.star Rotation matrix with gene sets in its columns and repetitions (rotations) in its rows. Note: ROws for ROtations.
#' @param left_sided Logical, defaults to TRUE. Are smaller p-values extremer?
#' @param filename A filename to do something. TODO what?
#' @return TODO please explain
#' @author Mathias Fuchs
do.MHC.R =  function(S.star, left_sided = TRUE, filename){     ## all functions for multiple hypothesis correction
 x = S.star  #TODO change this
 print(dim(x))
 orig = x[1, ]
 rot.mat = x[-1, ]

 tian = do.tian(orig.p = orig, rot.mat = rot.mat, gs.in.col = TRUE, left_sided = left_sided)
 MHC = list()
 MHC$tt = data.frame(
 #		geneset = colnames(x),
 orig = orig,
 marg.p.value =   tian$marg.p.value,
 q.value = 	 tian$q.value,
 q.least_mt_maj = tian$lmm,
 q.gr_mt_maj    = tian$gmm,
 ExpNrFD =        tian$ExpNrFD,
 NrDiscoveries =  tian$NrMoreExtreme,
 BH.q.value = p.adjust(tian$marg.p.value, method = "BH"),
 BY.q.value = p.adjust(tian$marg.p.value, method = "BY")
 );
 MHC$ttu = MHC$tt
 MHC$tt = MHC$tt[order(tian$lmm,tian$NrMoreExtreme), ]
 MHC$pi0 = tian$pi0
 if (!missing(filename)){
  require(R2HTML)
  require(xtable)
  HTML(MHC$tt, paste("results_html/", filename, ".html", sep = ""))
  system(paste("chromium-browser", paste("results_html/", filename, ".html", sep = "")))
  print(xtable(MHC$tt), type = "latex", file = paste("results_tex/", filename, ".tex", sep = ""), digits = 6, append = FALSE, floating = FALSE)
  return()
 }
 else {
  return(MHC)
 }
}

#############################
### Handy small functions ###
#############################

#' Test an Expression Matrix with Limma. Works with any design.
#' @param X Expression matrix with samples in columns and genes in rows.
#' @param group Response vector.
#' @return Fit object from `limma' function `eBayes'.
#' @author Stephan Artmann
limma.test = function (X,group) {
  design = model.matrix(~group);
  fit = lmFit(X, design)
  fit = eBayes(fit)
  fit;
}


#' Same as limma.test, but returns only the moderated t-statistics.
#' @param X Expression matrix with samples in columns and genes in rows.
#' @param group Response vector.
#' @return Vector of moderated t-statistics.
#' @author Mathias Fuchs
limma.test.slim = function(X, group) {
  fit = lmFit(X, model.matrix(~ group))
  fit$coefficients / fit$stdev.unscaled / sqrt(squeezeVar((fit$sigma)^2, fit$df.residual)$var.post)
}


#' Same as limma.test.slim, but returns ordinary instead of moderated t-statistics.
#' @param X Expression matrix with samples in columns and genes in rows.
#' @param group Response vector.
#' @return Vector of ordinary t-statistics.
#' @author Mathias Fuchs
limma.test.ordinaryt = function(X, group) {
  fit = lmFit(X, model.matrix(~ group))
  fit$coefficients / fit$stdev.unscaled / fit$sigma
}


# ' Simple implementation of the t-statistic assuming equal variances in both groups.
# ' Works so far only in two group case!
# ' @author Stephan Artmann
# t.test.simple = function(X,group) {
# r = factor(group);
# m = apply(X[,r == levels(r)[1]],1,mean) - apply(X[,r == levels(r)[1]],1,mean);  ### Differences in Mean
# v = apply((X - apply(X,1,mean))^2,1,sum)/(ncol(X)-1);                           ### Calculate Variance
# m/sqrt(v)/sqrt(1/length(which(r == levels(r)[1])) + 1/length(which(r == levels(r)[2])))
#}

#' Calculate one-sided p-values given two-sided ones and the fold change.
#' @param p Two-sided p-values.
#' @param M M-value, positive for up-, negative for down-regulated genes.
#' @param lower Logical, shall one-sided p-values stand for down-regulation?
#' @return One-sided p-values.
#' @author Stephan Artmann
p.one.sided = function(p,M,lower=FALSE) {
if (lower) M = -M;
p.star = p/2;
p.star[M < 0] = 1 - p.star[M<0];
p.star
}


#' Internal algorithm: Make limma test one-sided
#' @param fit Result of "lmFit" and "eBayes" functions in "limma" package.
#' @param lower Shall one-sided p-value indicate down-regultation?
#' @return One-sided p-values
limma.one.sided = function (fit,lower=FALSE) {
  se.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled
  df.total <- fit$df.prior + fit$df.residual
  pt(fit$t, df=df.total, lower.tail=lower)[,2]
}


#' Internal Function to test Gene Sets.
#' @param X Expression matrix.
#' @param response Response vector.
#' @param gene_sets_input Data.frame of gene sets.
#' @param nrot Number of rotations to be performed.
#' @param tests Which tests shall be applied? Currently, Wilcoxon ("W") or Goeman ("G") can be chosen. Put "all" if you want all tests to be performed.
#' @param one.sided Shall one-sided results be returned as well? So far only possible for Wilcoxon statistic.
#' @return List with matrix of rotated statistics and other fancy stuff.
#' @author Stephan Artmann
test.gene.sets = function (X, response, gene_sets_input, nrot = 1000 ,tests=c("W","G"),one.sided=TRUE, verbose=FALSE) {
 ### Check Input ###
 if(is.null(rownames(X))) stop ("Please specify row names of X");
 if (length(tests) == 0) stop("Please provide a gene set test: Either W for the Wilcoxon statistic, or G for the Goeman statistic.")
 if (all(tests == "all")) tests = c("W","G"); 
 if ("W" %in% tests & one.sided) tests = c(tests,"W.l","W.h");
 all.G = FALSE;
 if (all(tests == "G")) all.G = TRUE;
 testNo = length(tests);

 ### Load necessary libraries ###
 library(limma)


 if (length(levels(as.factor(response))) == 2) {  ### only do the following if the response has two levels:
  response = response - mean(response)             ### response has mean zero, and will keep it
  response = response / sqrt(sum(response^2))      ### likewise with standard deviation
 }
 gene_sets_input = as.data.frame(gene_sets_input)
 sizes = c()

 X = data.matrix(X)
 m = length(response)

 stopifnot(dim(X)[2] == m)

 geneSets = as.vector(unique(gene_sets_input[,2]))
 gsNo = length(geneSets)

 ### Original p-value ###
 S = list();                                               ### S is the list of statistics applied
 for (i in 1:testNo) {
  S[[i]] = c(NA);
 }
 index = list()
 others = list()

 ### the preparation
 modt = limma.test.slim(X,response);
 ranks = rank( -abs(modt), ties.method="random");
 if (one.sided) {
   fit = limma.test(X, response)
   ranks.l = rank(limma.one.sided(fit,lower=TRUE),ties.method="random");
  ranks.h = rank(-(ranks.l));
 }


 if (all.G) Lambda = list();
 mappable_indices = c();
 j1 = 0
 for (j2 in 1:gsNo) { ### CAUTION - why doesnt it work for just a few j_2 ?
   print(geneSets[j2])
   genes = as.vector(gene_sets_input[which(gene_sets_input[ ,2] == geneSets[j2]), 1]);   # here, "which" was missing.
   indx = match(as.character(genes), rownames(X));
   indx = indx[!is.na(indx)]	
  if (length(indx) > 0){
   j1 = j1 + 1
   mappable_indices = c(mappable_indices, j2)
   othrs = rep(TRUE, length(ranks))
   othrs[indx] = FALSE
   othrs = which(othrs)  
   index[[j1]] = indx
   others[[j1]] = othrs
   if (all.G) {
    if (length(response) < length(index[[j1]])) {
     Lambda[[j1]] = svd(X[index[[j1]], ])$d;
     S [[1]] [j1] = sum((Lambda[[j1]] * response)^2);
    } else {
     print(TRUE)
     Lambda[[j1]] = NULL; ##TODO hier noch ueberlegen, wie man diesen Fall abfangen will !!!s
    }
   } else {
    if ("W" %in% tests) {
     S[[which(tests == "W")]] [j1] = calculate.U(ranks, indx); 
    }
    if ("W.l" %in% tests) {
     S[[which(tests == "W.l")]] [j1] = calculate.U(ranks.l, indx); 
    }
    if ("W.h" %in% tests) {
     S[[which(tests == "W.h")]] [j1] = calculate.U(ranks.l, indx); 
    }
    if ("G" %in% tests) {
     S[[which(tests == "G")]] [j1] = Q(X[indx,], response);           ### TODO schau bitte nochmal drueber, ob das Deinem Aufruf entspricht (stephan)
                                                                      ### TODO sorge bitte dafuer, dass bei Q smaller.is.extremer=TRUE rauskommt (im Zweifel mit -1 multiplizieren - danke! (stephan)
    }
   }
   sizes = c(sizes, length(indx))
  }
 }

 geneSets = geneSets[mappable_indices]
 gsNo = length(mappable_indices)

 S.star = list();
 for (i in 1:testNo) {
  S.star[[i]] = rep(NA, nrot*gsNo);	
  dim(S.star[[i]]) = c(nrot,gsNo);
 }

 ### Rotation ###
 for (i in 1:nrot) { ## core engine
  if (verbose) print(paste("rotation number", i))
  if (all.G) {
	  rr = random_response(m); 
  } else {
   X.star = X %*% rotation(length(response), method = 2)   
   modt = limma.test.slim(X.star,response)
   ranks = rank(-abs(modt), ties.method = "random")
   if (one.sided) {
    fit = limma.test(X.star, response)
    ranks.l = rank(limma.one.sided(fit,lower=TRUE),ties.method="random");
    ranks.h = rank(-(ranks.l));
   }
  }
  for (j in 1:gsNo) {
   if (all.G) {
    if (is.null(Lambda[[j]])) next;
    S.star[[1]][i, j] = sum((Lambda[[j]]*rr)^2);
   } else {
    if ("W" %in% tests) S.star[[which(tests == "W")]][i, j] = calculate.U(ranks, index[[j]]);
    if ("W.l" %in% tests) S.star[[which(tests == "W.l")]][i, j] = calculate.U(ranks.l, index[[j]]);
    if ("W.h" %in% tests) S.star[[which(tests == "W.h")]][i, j] = calculate.U(ranks.h, index[[j]]);
    if ("G" %in% tests) S.star[[which(tests == "G")]][i,j] = Q(X.star[index[[j]],],response);
   }
  }
 }
 for (i in 1:length(S.star)) {
  S.star[[i]] = rbind(S[[i]],S.star[[i]]);
  colnames(S.star[[i]]) = geneSets;
 }
 names (S.star) = tests;
 x = list()
 x$S.star = S.star
 x$sizes = sizes
 return(x)
}

###############################
### Current Test Statistics ###
###############################

#' Calculate Mhann-Whitney's U.
#' @param ranks Ranks of p-values.
#' @param index Logical or numeric vector indicating where the gene set is.
#' @return U1 statistic of Mhann-Whitney.
#' @author Stephan Artmann
calculate.U = function(ranks,index) {
 n1 = length(index);
 n2 = length(ranks) - n1;
 U1 = sum(ranks[index]) - n1*(n1+1)/2;
# U2 = n1*n2-U1;
 U1 = U1 - n1 * n2 / 2; # subtract mean
 U1 = U1 / sqrt(n1*n2*(n1+n2 +1)/12) # divide by standard deviation
return(U1);

}


#' Reimplementation of Goeman's Q.
#' @param X Expression matrix.
#' @param response Response vector.
#' @return Goeman's Q.
#' @author Mathias Fuchs
Q = function (X, response){
    response = response - mean(response)
    response = response / sqrt(sum(response^2))

    stopifnot (length(response) == dim(X)[2])
   for (i in 1:dim(X)[1]){
	X[i, ] = X[i, ] - mean(X[i, ])
#	X[i, ] = X[i, ] / sqrt(sum(X[i, ]^2))
    }
    return(as.numeric(t(response) %*% t(X)%*% X %*% response) / dim(X)[1])
}
##########################
### Rotation Functions ###
##########################

#' Gives a vector a unique spherical vector in the subspace of zero mean, of length 1
#' @param m Probably length? TODO
#' @author Mathias Fuchs
random_response = function(m) { 
# manufacture a unit vector in dimension m-1
    stopifnot(m >= 2)
    x = rnorm(m)
    x = x - mean(x)
    x = x / sqrt(sum(x^2))
}

#' Generates a random rotation matrix according to the nested subgroup method described Diaconis Shashahani 1987
#' @param n Number of replicates
#' @param method method = 1 generates full rotations, method = 2 generates a rotation in the subgroup of mean-preserving rotations.
#' @param genes_in_rows Logical, are genes in rows?
#' @author Mathias Fuchs
rotation = function (n, method, genes_in_rows ) {
    stopifnot (n >= 2)
    if (method == 1){
	if (n == 2){
	    angle = runif(1) * 2 * pi
	    s = sin(angle)
	    c = cos(angle)
	    A = matrix(NA, 2, 2)
	    A[1, 1] = s
	    A[2, 1] = c
	    A[1, 2] = -c
	    A[2, 2] = s
	    if(runif(1) < .5){    # the determinant also has to be random
		    A = A %*% diag(c(-1, 1))
	    }
	    return(A)
	}
	else {
	    gamma = c() # we use the algorithm of Diaconis  Shashahani 1987, page 23 - yeah 8)
	    gamma = rotation(n - 1, method = 1) # this is kind of recursive! - no, it IS recursive :)
	    gamma = cbind(rep(0, n-1), gamma)
	    gamma = rbind(c(1, rep(0, n-1)), gamma)
	    randomspherepoint = rnorm(n)
	    randomspherepoint = randomspherepoint / sqrt(sum(randomspherepoint^2))
	    eoneminusv = c(1, rep(0, n-1)) - randomspherepoint 
	    x = eoneminusv / sqrt(sum((eoneminusv)^2))
	    dim(x) = c(n, 1)
	    return((diag(n) - 2 * x %*% t(x)) %*% gamma)
	}
    }
    else {
	if (n == 2){
	    stop("method 2 makes no sense in dimension 2")
	}
	else {
	    turn = matrix(NA, nc = n, nr = n)
	    M = c()
	    M = rep(1, n) / sqrt(n) # this vector is orthogonal to what we need
	    dim(M) = c(n, 1)
	    turn = c() # this matrix rotates in such a way that the last vector becomes M
			# the explicit formula used here, is quite funny.
			# it arises from solving some quadratic equation on the entries.
		turn = matrix(rep((sqrt(n) - 1)/(n-1), n^2), n, n)
		turn[,1] = rep(1, n)
		turn[1,] = rep(1, n)
		diag(turn) = c(1, rep( -(1 + (sqrt(n) - 1) * (n-2) / (n-1)) , n-1))
		turn = turn %*% diag(rep(1 / sqrt(n), n))
	    r = c()
	    r = rotation(n-1, method = 1) # this only seems to be recursive, but isn't, because it's the other method.
	    r = cbind(rep(0, n-1), r)
	    r = rbind(c(1, rep(0, n-1)), r)
	    return(turn %*% r %*% t(turn))
	# we first map the subspace of zero mean to the space othogonal to the first base vector, then we rotate the last n-1 coordinates, then we turn back
	}
    }
}

#############################################
### Functions to deal with .star matrices ###
#############################################


#' Calculation of the least monotone majorant of q-values.  
#' @param q Vector of q-values.
#' @param p Vector of p-values. Necessary to monotonise the q-values according to p-values.
#' @return Least monotone majorant as vector.
#' @author Stephan Artmann, Mathias Fuchs
least_monotone_majorant = function(p, q){
  stopifnot(is.vector(q))
  stopifnot(is.numeric(q))
  stopifnot(length(p) == length(q))

  lq = length(q)

  i <- lq:1
  o <- order(p,1-q, decreasing = FALSE)
  ro <- order(o)
  lmm = pmin(1, cummax(q[o]))[ro]

  return(lmm)
}

#' Calculation of the greatest monotone minorant of q-values. TODO change name
#' @param q Vector of q-values.
#' @param p Vector of p-values. Necessary to monotonise the q-values according to p-values.
#' @return Greatest monotone minorant as vector.
#' @author Stephan Artmann, Mathias Fuchs
gr_monotone_majorant = function(p, q){
  stopifnot(is.vector(q))
  stopifnot(is.numeric(q))
  stopifnot(length(p) == length(q))

  lq = length(q)

  i <- lq:1
  o <- order(p,1-q, decreasing = TRUE)
  ro <- order(o)
  lmm = pmin(1, cummin(q[o]))[ro]

  return(lmm)
}

#' Tian mutliple testing correction. TODO remove loop
#' @param orig.p Original p-value vector.
#' @param rot.mat Matrix of rotated gene set statistics.
#' @param gs.in.col logical, where are the gene sets in rot.mat?
#' @param left_sided logical, lead lower statistics to lower p-values?
#' @return TIAN, what else? TODO
#' @author Stephan Artmann, Mathias Fuchs
do.tian = function(orig.p,rot.mat, gs.in.col=TRUE, left_sided = TRUE) {
	stopifnot(is.numeric(rot.mat))
	stopifnot(is.numeric(orig.p))

 rot.mat.rank = function (orig.p,rot.mat) {
  rank.orig = rank(orig.p,ties.method="max");
  dim(rot.mat) = NULL;
  M = rank(c(orig.p,rot.mat),ties.method="max");
  M[1:length(orig.p)] - rank.orig;
 }


    require(qvalue)
    if (!gs.in.col) rot.mat = t(rot.mat);
	stopifnot(length(orig.p) == dim(rot.mat)[2])

	if(!left_sided) {
	 rot.mat = (-1) * rot.mat
	 orig.p = (-1) * orig.p
	}

    marg.p = c()
    for (i in 1:length(orig.p)){
	marg.p = c(marg.p, marginal_p(orig.p[i], rot.mat[, i], left_sided = left_sided)) ## the marginal p-values
    }
    pi0 = try(qvalue(marg.p)$pi0)
    if (inherits(pi0, "try-error")) {pi0 = 1}
    print(paste("Estimated portion of true nulls =", pi0))

    #NrMoreExtreme = rank(marg.p, ties = "max")       ## this is the number of more extreme p-values, i.e. the number of discoveries
    NrMoreExtreme = rank(orig.p, ties = "max")       ## this is the number of more extreme p-values, i.e. the number of discoveries
    ExpNrFD = c()

    #for (i in 1:length(orig.p)){
    #  rot.mat[, i] = rank(rot.mat[, i])
    #}			




        ExpNrFD = rep(NA, length(orig.p))
        ExpNrFD = pi0 * rot.mat.rank(orig.p,rot.mat)/nrow(rot.mat)

  ## this is the number of false discoveries i. e. the portion of true nulls, times the null-expected number of more extreme test statistics
    
    TIAN = list()

    TIAN$marg.p.value = marg.p
    TIAN$q.value = ExpNrFD / NrMoreExtreme   ## the q-value is the number of false discoveries divided by the number of discoveries
    TIAN$pi0 = pi0
    TIAN$ExpNrFD = ExpNrFD
    TIAN$NrMoreExtreme = NrMoreExtreme 	
    TIAN$lmm = least_monotone_majorant(TIAN$marg.p.value, TIAN$q.value)
    TIAN$gmm = gr_monotone_majorant(TIAN$marg.p.value, TIAN$q.value)
    return(TIAN)
}


#' Generate marginal p-values.
#' @param x a real number interpreted as the value of a test statistic
#' @param y a vector interpreted as an empirical sample distribution
#' @param left_sided the notorious smaller-is-extremer-logical
#' @author Stephan Artmann, Mathias Fuchs
marginal_p = function (x, y, left_sided = TRUE){
    if (!left_sided) {
	x = -x
	y = -y
    }
    return ((length(which(y <= x)) + 1) / (length(y) + 1))
}



###############################################################
### Additional Functions ###
###############################################################

#' do.MHC.file
#' @param filename filename
#' @param left_sided left_sided
#' @return return
#' @author Mathias Fuchs
do.MHC.file =  function(filename, left_sided = TRUE){ ## all functions for multiple hypothesis correction
    x = data.matrix(read.csv(filename))
    return(do.MHC.R(x, left_sided = left_sided))
}

#' Data an auxiliary function designed to load data from a subdirectory or from GEO
#' @param fn optional, a downloaded series matrix filename, possibly with path
#' @param GSEid optional, a GSE identifier
#' @param two_class logical, is it only a two-class-problem?
#' @return returns list with expression and, case given, differential information
#' @author Mathias Fuchs
Data = function(fn = NA, GSEid = NA, two_class = FALSE) {
  require(GEOquery)
  if (file.exists(fn))
    xx = getGEO(file = fn)
  else
    xx = getGEO(GSEid)
  x = list()
  x$m  = exprs(xx)
  x$m = Strip_NA_rows(x$m)
  if (two_class){
      require(limma)
    x$DA.m = limma.test(x$m, x$group)
    print("significant p-values")
    print(length(which(x$DA.m$p[,2] < 0.05)))
    print("out of")
    print(length(x$DA.m$p[,2]))
    print("significant mRNA q-values")
    q = p.adjust(x2$DA.m$p[,2], method = "BH")
    print(length(which(q < 0.05)))
    x$m.short = x$m[order(x$DA.m$p[,2])[1:150],]
  }
  return(x)
}

#' Strip_NA_rows
#' @param A a data matrix
#' @return returns the same matrix without those rows with at least one NA value
#' @author Mathias Fuchs
Strip_NA_rows = function(A){
   stopifnot(is.matrix(A))
   A = data.matrix(A)
   index = c()
   for (i in 1:dim(A)[1]) {
      if (any(is.na(A[i,]))) {
	  index = c(index, i)
      }
    }
    if (length(index) > 0){
      A = A[- index, ]
    }
    return(A) 
}

#' visualize TODO check spelling, visualize -> visualise :P
#' @param A input expression matrix
#' @param response the usual group/class vector
#' @param method one out of "scatterplot3d", "rotate", "pairs" 
#' @param NoPCA_for_pairs the number of principal components for the pairs method
#' @return returns nothing, just tries to be beautiful
#' @author Mathias Fuchs
visualize = function(A, response, method = "rotate", NoPCA_for_pairs= 3) {
    stopifnot(is.matrix(A))
    stopifnot(length(response) == dim(A)[2])
    responses = unique(response)
    palette(rainbow(length(responses))) # we need less colors, but this produces nicer ones
    cols = rep(NULL, length(response))
    for (i in 1:length(responses)){
	cols[which(response == responses[i])] = palette()[i]
    }
    if (method == "scatterplot3d"){
	require(scatterplot3d)
	scatterplot3d(t(ShortenMatrixOntoPrincipalComponent(A,3)), color = cols)
    }
    else if (method == "rotate"){
	require(rgl)
	plot3d(t(ShortenMatrixOntoPrincipalComponent(A,3)), col = cols, type = "s", , aspect = FALSE)
    }
    else if (method == "pairs"){
	pairs(t(ShortenMatrixOntoPrincipalComponent(A, NoPCA_for_pairs)), col = cols)
    }
}

#' Function.
#' @param x x
#' @param n_new n_new
#' @author Mathias Fuchs
ShortenMatrixOntoPrincipalComponent = function(x, n_new){
  stopifnot(is.matrix(x))
  stopifnot(n_new > 0)
  n_old = dim(x)[1]
  m = dim(x)[2]
  if (n_new == n_old) return(x) # then there is nothing to do, anyway, this should not happen
  if (n_new > n_old) stop ("There is nothing to shorten")
  svd_x = svd(x)
  smallest_sv = min(svd_x$d[which(svd_x$d > 0)])
  svd_x$d[which(svd_x$d == 0)] = smallest_sv    # in order to prevent numerical accidents where a singular value is zero
    if (n_new == 1){
	return(svd_x$d[1] * t(svd_x$v[,1]))
    }
    else if (n_new < m){                           # this is the typical case
	return(diag(svd_x$d[1:n_new]) %*% t(svd_x$v[,1:n_new]))
    }
    else if (n_new == m){
	reurn(diag(svd_x$d) %*% t(svd_x$v))
    }
    else if (n_new > m){
	stop("Wrong dimensions in svd.")
    }
}

# end of file

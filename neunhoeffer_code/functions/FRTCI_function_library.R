##
##  A "Package" for doing Permutation Test inference for heterogeneous
##  treatment effects.
##
## To use, source this file at the top of your analysis file.  See
## the demo analysis script for a specfic example.
##
## (C) 2015.   Peng Ding, Avi Feller, Luke Miratrix




################################################################
###
### A library of test statistics
###
### All are a function of Y, Z for first two arguments
###
################################################################


# Calculate classic (not shifted) KS statistic
# Code modified version of R's ks.test 
# If tau passed, Y1 will be shifted by tau.
KS.stat <- function( Y, Z, tau = NULL, alternative = c("two.sided", "less", "greater") ) {
    x = Y[Z==1]
    y = Y[Z==0]
    if ( !is.null( tau ) ) {
        x = x - tau
    }
    alternative <- match.arg(alternative)
    
    n.x <- length(x)
    n.y <- length(y)
    
    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
    if (length(unique(w)) < (n.x + n.y)) {
        z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
    }
    STATISTIC <- switch(alternative, two.sided = max(abs(z)), 
                        greater = max(z), less = -min(z))
    
    STATISTIC
}


## Shifted kolmogorov-smirnov statistic
## Calculate KS distance between Y0 and Y1 shifted by estimated tau.
SKS.stat <- function(Y, Z)
{
    Y1 = Y[Z==1]
    Y0 = Y[Z==0]
    
    Y1.star   = Y1 - mean(Y1)
    Y0.star   = Y0 - mean(Y0)
    
    unique.points = c(Y1.star, Y0.star)
    
    Fn1 = ecdf(Y1.star)
    Fn0 = ecdf(Y0.star)
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
    
}


## Shifted kolmogorov-smirnov statistic with covariates
## to increase precision.
##
## This is the test statistic used in the JRSS B Paper.
SKS.stat.cov.pool <- function(Y, Z, X)
{
    this.lm <- lm(Y ~ Z + X)
    Y1.star <- this.lm$res[Z == 1]
    Y0.star <- this.lm$res[Z == 0]
    
    unique.points = c(Y1.star, Y0.star)
    
    Fn1 = ecdf(Y1.star)
    Fn0 = ecdf(Y0.star)
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
}


## Shifted kolmogorov-smirnov statistic with covariates
## with model for outcomes calculated on control group only.
## This avoids "splitting" the treatment variation between tx 
## and co groups.  
## We recommend this method over the "pool" method.
SKS.stat.cov <- function(Y, Z, X)
{
        this.lm <- lm(Y ~ X, subset = Z == 0)
        
        Y0.hat <- predict(this.lm, newdata = as.data.frame(X))
        Y0.res <- Y - Y0.hat
        
        Y1.star <- Y0.res[Z == 1] - mean(Y0.res[Z == 1])
        Y0.star <- Y0.res[Z == 0] - mean(Y0.res[Z == 0])
        
        unique.points = c(Y1.star, Y0.star)
        
        Fn1 = ecdf(Y1.star)
        Fn0 = ecdf(Y0.star)
        
        difference = Fn1(unique.points) - Fn0(unique.points)
        
        return(max(abs(difference)))
}





## Shifted kolmogorov-smirnov statistic with a linear treatment
## effect model defined by W
##
## This will attempt to remove any systematic variation corresponding
## to W and then return a SKS statistic on the residuals to measure
## any variation "left over".
##
## X are _additional_ covariates to adjust for beyond those involved
## in treatment effect model.  It will automatically ajust for W as
## well.  Do not put a covariate in for both X and W.
##
## For use in the FRTCI.interact method
## This is the test statistic used in the JRSS B Paper.
SKS.stat.int.cov.pool <- function( Y, Z, W, X=NULL )
{ 
    if ( !is.null( X ) ) {
        this.lm <- lm( Y ~ Z + X + W + Z:W )
    } else {
        this.lm <- lm( Y ~ Z + W + Z:W )        
    }
    Y1.star <- this.lm$res[Z == 1]
    Y0.star <- this.lm$res[Z == 0]
    
    unique.points = c(Y1.star, Y0.star)
    
    Fn1 = ecdf(Y1.star)
    Fn0 = ecdf(Y0.star)
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
}



## Shifted kolmogorov-smirnov statistic with a linear treatment
## effect model defined by W
##
## This will attempt to remove any systematic variation corresponding
## to W and then return a SKS statistic on the residuals to measure
## any variation "left over".
##
## X are _additional_ covariates to adjust for beyond those involved
## in treatment effect model.  It will automatically ajust for W as
## well.  Do not put a covariate in for both X and W.
##
## For use in the FRTCI.interact method
## This is the test statistic used in the JRSS B Paper.
##
## This method first adjusts for baseline and then models treatment effect
## on the residuals to not split treatment effects.
##
## For use in the FRTCI.interact method.
## We recommend this method over the "pool" method.
SKS.stat.int.cov <- function( Y, Z, W, X=NULL )
{ 
    # First wipe out Y0 predicted by X via linear model
    if ( !is.null( X ) ) {
        this.lm <- lm(Y ~ X, subset = Z == 0)
        
        Y0.hat <- predict(this.lm, newdata = as.data.frame(X))
        Y <- Y - Y0.hat
    }        

    # now model treatment effect
    this.lm <- lm( Y ~ Z + W + Z:W )        

    Y1.star <- this.lm$res[Z == 1]
    Y0.star <- this.lm$res[Z == 0]
    
    unique.points = c(Y1.star, Y0.star)
    
    Fn1 = ecdf(Y1.star)
    Fn0 = ecdf(Y0.star)
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
}



##
## Other possible test statistics
##


## Shifted kolmogorov-smirnov statistic with covariates and rq
SKS.stat.cov.rq <- function(Y, Z, X)
{
    require( quantreg )
    this.lm <- lm(Y ~ Z + X)
    this.rq <- rq(Y ~ Z + X, tau = seq(0.05, 0.95, by = 0.05))
    
    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )
    
}


## Kolmogorov-smirnov statistic via quantile regression with covariates
rq.stat <- function(Y, Z, rq.pts = seq(0.1, 0.9, by = 0.1))
{
    require( quantreg )
    
    this.lm <- lm(Y ~ Z)
    this.rq <- rq(Y ~ Z, tau = rq.pts)
    
    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )
    
}


## Kolmogorov-smirnov statistic via quantile regression with covariates
## Conditional approach; see Koenker and Xiao (2002)
rq.stat.cond.cov <- function(Y, Z, X, rq.pts = seq(0.1, 0.9, by = 0.1))
{
    require( quantreg )
    
    this.lm <- lm(Y ~ Z + X)
    this.rq <- rq(Y ~ Z + X, tau = rq.pts)
    
    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )
    
}



## Kolmogorov-smirnov statistic via quantile regression with covariates
## Unconditional approach; see Firpo (2007)
rq.stat.uncond.cov <- function(Y, Z, X, rq.pts = seq(0.1, 0.9, by = 0.1))
{
    require( quantreg )
    
    ## propensity score model
    this.glm <- glm(Z ~ X, family = binomial(link = logit))
    
    pscore <- predict(this.glm, type = "response")
    ipw.weights <- ifelse(Z, 1/pscore, 1/(1 - pscore))
    
    
    this.lm <- lm(Y ~ Z, w = ipw.weights)
    this.rq <- rq(Y ~ Z, tau = rq.pts, w = ipw.weights)
    
    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )
    
}


##
## Test statistics for categorical covariates as discussed in the
## JRSS B Paper.
##

## Weighted average of the group-level SKS statistics
## @param W is a factor or categorical covariate of some sort.
t.WSKS = function( Y, Z, W ) {
    
    dd = ddply( data.frame(Y=Y,Z=Z,W=W), "W", summarize, 
                t.sks = SKS.stat( Y, Z ),
                n.k = length(Y) )
    n = length(Y)
    sum( dd$t.sks * dd$n.k / n )
}


## Subtract off group level treatment effect estimates and then look
## at KS statistic on residuals.  
##
## Distinct from the interacted lm in that the control units are not
## shifted and centered with respect to eachother.
t.SKS.pool = function( Y, Z, W ) {    
    dat <- data.frame( Y=Y, Z=Z, W=W )
    mns = ddply( dat, .(W, Z), summarize, mean=mean(Y) )
    taus = with( mns, mean[Z==1] - mean[Z==0] )
    KS.stat( dat$Y - dat$Z*taus[dat$W], dat$Z )
}





################################################################
##
## Sampling from confidence region of multivariate normal
##
################################################################

# Generate truncated multivariate normal draws
# where we don't return anything in the tails defined by alpha.
rcrmtvnorm <- function(n, mu, sigma, alpha)  {
    ## dimension
    K = length(mu)
    
    ## return a matrix
    Mat = matrix(0, n, K)
    
    ## critical value of chisq distribution of dof = K
    ctv.chisq = qchisq(p = 1 - alpha, df = K)
    
    count = 0
    repeat{
        
        draw.try = rmvnorm(n = 1, mean = mu, sigma = sigma)
        draw.try = as.vector(draw.try)
        
        ##within confidence region or not?
        chisq.try = t(draw.try - mu)%*%solve(sigma, draw.try - mu)
        if(chisq.try <= ctv.chisq)
        {
            count = count + 1
            Mat[count, ] = draw.try
            
            if(count == n) break       	
        }## if      	
        
    }## repeat	
    
    return(Mat)
}





################################################################
##
##  The conducting-a-test code
##
################################################################


# Generate a sequence of tau values in the confidence interval of 
# effects to search over with the permutation test.
get.tau.vector = function( Y, Z, X = NULL, gamma=0.001, grid.size=21, grid.gamma = 100*gamma ) {

    if ( is.null( X ) ) {
        
        te.hat <- mean(Y[Z == 1]) - mean(Y[Z == 0])
        te.se <- sqrt( var(Y[Z == 1])/sum(Z) + var(Y[Z == 0])/sum(1 - Z))
    
    } else {
        lm.tau <- lm( Y ~ Z + X )
    
        te.hat <- as.numeric(coef(lm.tau)["Z"])
        te.se <- as.numeric(sqrt( diag(vcov(lm.tau))["Z"]))
    }

    # oversample points near the estimated tau-hat
    te.MOE <- qnorm(1 - gamma/2)*te.se
    te.vec <- te.hat + (te.MOE/qnorm(1-grid.gamma)) * qnorm( seq( grid.gamma, 1-grid.gamma, length.out=grid.size ) )
    te.vec
        
    attr( te.vec, "te.hat" ) <- te.hat
    attr( te.vec, "te.se" ) <- te.se
    attr( te.vec, "te.MOE" ) <- te.MOE
    
    te.vec
}




# Calculate a set of different models for effects based on a linear model of W, controlling
# for X (and W).  These models all correspond to estimated nusiance parameters beta.
# Each possible beta gives a collection of different models of effects.
#
# I.e., estimate the coefficients a, d for
# Y ~ a Z + b X + c W + d W:Z
#
# Then for each model, calculate individual imputed treatment effects for all observations and the
# associated science tables of imputed potential outcomes.
#
# @param  Y, Z, X, W: outcome, treatment assignment, covariates, and treatment varying covariates
#                  X and W can be the same for a fully interacted linear model.
# @return
#        te.grid    grid.size x p sized matrix with each row j corresponding to a different model of
#                   effects with specific beta.  p is number of columns in W.
#        Y0.mat, Y1.mat   Two N x grid.size matrices.  Each column j corresponds to a specific model of 
#                   effects with specific beta value.  Each column j corresponds to the same row in te.grid
#         te.mat    N x grid.size matrix with each row i corresponding to treatment effect for 
#                   unit i under model of effects j.
get.testing.grid = function( Y, Z, W, X=NULL, gamma=0.0001, grid.size=150 ) {
    
    ## get sample of treatment effect models to calculate p-values for
    if ( !is.null(X) ) {
        lm.tau <- lm(Y ~ Z + W + Z:W + X)
        # if have duplicates with X and W, need to drop from consideration
        drp = is.na( coef(lm.tau) )
        cof <- coef(lm.tau)[ !drp ]                
    } else {
        lm.tau <- lm( Y ~ Z + W + Z:W )
        cof <- coef(lm.tau )
    }
    te.model = grep( "Z", names(cof) )
        
    stopifnot( all( rownames( vcov(lm.tau) ) == names(cof) ) )
    
    te.hat <- as.numeric( cof[te.model] )    
    te.cov <- vcov(lm.tau)[te.model,te.model]
    
    te.grid = te.hat
    if ( grid.size > 1 ) {
        # sample points from our confidence region, focusing on points
        # close to our point estimate
        te.grid <-  rbind( te.grid, rcrmtvnorm(grid.size - 1, mu = te.hat, sigma = te.cov, alpha = gamma ) )
    }
    
    # calculate individual treatment effects for different treatment impact models
    # each row is individual and each column is specific treatment effect model
    te.mat <- t(te.grid %*% t(cbind(1, W))) 
    
    Y0.mat <- Y*(1 - Z) + (Y - te.mat) * Z
    Y1.mat <- Y*Z + (Y + te.mat) * (1 - Z)
        
    list( te.grid=te.grid, te.mat=te.mat, Y0.mat=Y0.mat, Y1.mat=Y1.mat )
}







# Utility function used by FRTCI and FRTCI.interact to actually generate the 
# permutations.  Don't use directly.
generate.permutations = function( Y, Z, test.stat, Y0.mat, Y1.mat, B, get.z.star=NULL, verbose = TRUE, ... ) {
    
    ## SET UP STORAGE MATRICES
    n.te.vec <- ncol(Y0.mat)
    ks.mat <- matrix(NA, nrow = B, ncol = n.te.vec)
    
    ## CALCULATE OBSERVED TEST STATISTICS
    ks.obs <- test.stat( Y, Z, ... )
    
    if ( verbose ) {
        pb <- txtProgressBar(min = 0, max = B, style = 3)
    }
    
    ## RANDOMIZATION DISTRIBUTION
    for(b in 1:B){
        if ( verbose ) {
            setTxtProgressBar(pb, b)
        }
        if ( !is.null( get.z.star ) ) {
            Z.star = get.z.star( Z, ... )
        } else {
            Z.star <- sample(Z)
        }
        
        ##### CI method
        ci.out <- sapply(1:n.te.vec, 
                         function(i){ 
                             ## CALCULATE RANDOMIZED VALUE
                             Yobs.star <- Z.star*Y1.mat[,i] + (1 - Z.star)*Y0.mat[,i]
                             
                             ## COMPUTE TEST STATISTICS
                             ks.star <- test.stat(Yobs.star, Z.star, ... )
                             
                             ## OUTPUT
                             ks.star
                         })  
        
        ks.mat[b,] <- as.numeric(ci.out)
    }
    if ( verbose ) {
        close(pb)
    }
    
    ## CALCULATE P-VALUES
    ci.p <- apply(ks.mat, 2, function(ks.star){
        sum(ks.star >= ks.obs)/B
    })
    
    list( ks.obs=ks.obs, ks.mat=ks.mat, ci.p=ci.p )
}



## Conduct the FRT CI Method on passed data
##
## @param Y  Observed outcome vector
## @param Z  Treatment assigment vector (1=Tx, 0=Co)
## @param test.stat  Test statistic function to use on the data.
## @param B  Number of permutations to take.
## @param gamma How wide of a CI to make around tau-hat for search
## @param grid.gamma Parameter to govern where the grid points are sampled.  Bigger values
##    means more samples towards the estimated tau-hat.
## @param grid.size Number of points in the grid.
## @param te.vec Vector of taus to examine if you want to override generating ones automatically.
## @param return.matrix  TRUE means give back the matrix of all the imputed statistics.  FALSE do not.
## @param verbose  TRUE means print out progress bar when fitting and other diagnostics.
## @param ... Extra arguments passed to the generate.permutations function (and the test.stat function).
FRTCI <- function( Y, Z, test.stat = SKS.stat, B=500, 
                   gamma=0.0001, grid.gamma=100*gamma, 
                   grid.size=151,
                   te.vec=NULL, return.matrix=FALSE, 
                   verbose=TRUE, ... ) {
    
    
    if ( is.null(te.vec) ) {
        if ( grid.size %% 2 == 0 ) {
            grid.size <- grid.size+1
        }
        
        te.vec <- get.tau.vector( Y, Z, gamma=gamma, grid.size=grid.size, grid.gamma=grid.gamma )
    } else {
        grid.size = length( te.vec )
    }
    te.hat <- attr( te.vec, "te.hat" )
    te.se <- attr( te.vec, "te.se" )
    te.MOE <- attr( te.vec, "te.MOE" )
    
    ## IMPUTE MISSING POTENTIAL OUTCOMES
    Y1.mat <- sapply(te.vec, function(te) ifelse(Z, Y, Y + te) )
    Y0.mat <- sapply(te.vec, function(te) ifelse(Z, Y - te, Y) )
    
    res <- generate.permutations( Y, Z=Z, test.stat=test.stat, Y0.mat=Y0.mat, Y1.mat=Y1.mat, B=B, verbose=verbose, ... )
    
    ci.p = res$ci.p + gamma
    
    t = res$ks.obs
    p.value = max( ci.p )
    n=length(Y)
    method = paste( "FRT CI Test for Tx Effect Heterogeniety with ", substitute(test.stat), sep="" )
    DAT = paste( n, " observations", sep="")
    
    if ( !return.matrix ) {
        ks.mat=NULL
    } else {
        ks.mat = res$ks.mat
    }
    
    structure(list(statistic = t, p.value = p.value,
                   p.value.plug = ci.p[(grid.size+1)/2],
                   method=method,
                   data.name = DAT,
                   Y=Y, Z=Z, n=n, ci.p=ci.p, te.vec=te.vec,
                   te.hat=te.hat,
                   te.SE=te.se,
                   te.MOE = te.MOE,
                   B=B, gamma=gamma, ks.mat=ks.mat ),
              
              class = "FRTCI.test")
} 



# Return the plug-in p-value for plugging in tau-hat under the permutation method
# Just calls FRTCI.test with the specific, single tau of tau-hat.
FRTplug <- function( Y, Z, test.stat=SKS.stat, tau.hat=mean(Y[Z == 1]) - mean(Y[Z == 0]), ... ){
    mth = FRTCI( Y, Z, test.stat, te.vec=c(tau.hat), ...)
    mth$method = paste( "FRT Plug-in Test for Tx Effect Heterogeniety with ", substitute(test.stat), sep="" )
    mth
}







## Conduct the FRT CI Method, adjusting for covariates using a linear model, on passed data
## for a model of effects W'beta with unknown beta.
##
## See FRTCI() for further description of parameters.
##
## @param gamma How wide of a CI to make around tau-hat for search
## @param grid.gamma Parameter to govern where the grid points are sampled.  Bigger values
##    means more samples towards the estimated tau-hat.
## @param grid.size Number of samples for grid.
## @param test.stat  Test function.  Must take Y, Z, X, and W in that order.
## @param ... Extra arguments passed to the generate.permutations and test.stat function.
FRTCI.interact <- function( Y, Z, W, X=NULL, test.stat = SKS.stat.int.cov, B=500, 
                            gamma=0.0001, grid.gamma=100*gamma, 
                            grid.size=151, return.matrix=FALSE, 
                            verbose=TRUE, ... ) {
    
    grid.info = get.testing.grid( Y, Z, W=W, X=X, gamma, grid.size )
    
    te.MOE = NA
    
    te.grid = grid.info$te.grid
    Y1.mat <- grid.info$Y1.mat
    Y0.mat <- grid.info$Y0.mat
    
    res <- generate.permutations( Y, Z, test.stat, Y0.mat, Y1.mat, B=B, verbose=verbose, X=X, W=W, ... )
    
    t = res$ks.obs
    ci.p = res$ci.p + gamma
    p.value = max( ci.p )
    n=length(Y)
    method = paste( "FRT CI Test for Tx Effect Heterogeniety Beyond a Systematic Model with ", substitute(test.stat), sep="" )
    DAT = paste( n, " observations", sep="")
    
    if ( !return.matrix ) {
        ks.mat=NULL
    } else {
        ks.mat = res$ks.mat
    }
    
    structure(list(statistic = t, p.value = p.value,
                   p.value.plug = ci.p[1],
                   method=method,
                   data.name = DAT,
                   Y=Y, Z=Z, n=n, ci.p=ci.p, te.grid=te.grid,
                   B=B, gamma=gamma, ks.mat=ks.mat,
                   W=W, X=X ),
              
              class = "FRTCI.test")
} 










################################################################
##
##  Utility functions for dealing with result of test
##  code.
##
################################################################




# Print out FRTCI.test object with various information in pretty form.
print.FRTCI.test = function( x, ... ) {
    cat( "\t", x$method, "\n" )
    cat( "data: ", x$data.name, "\n" )
    cat( "t =", x$statistic, ", p-value =", x$p.value, " (plug = ", x$p.value.plug, ")\n" )
    rng = range(x$ci.p)
    cat( "\tp-value range =", rng[1], "-", rng[2] , "\n" )
    
    if ( !is.null( x$te.vec ) ) {    
        rng = range( x$te.vec )
        cat( "CI range = ", rng[1], "-",rng[2], "\tMOE =", x$te.MOE, "\n" )
        cat( "\tNeyman Difference point estimate of average =", x$te.hat, " SE =", x$te.SE,  "\n" )
    }
    cat( "\t# tested points =", length( x$ci.p ), "\n" )
    cat( "\tB = ", x$B, "\tgamma = ", x$gamma, "\n" )
    
    if ( !is.null( x$W ) ) {
        cat( "\tCorrected for", paste( colnames(x$W), sep=", " ) )        
    }

    if ( !is.null( x$X ) ) {
        cat( "\tAdjusted for", paste( colnames(x$X), sep=", " ) )        
    }
    
    
}





# Give confidence bounds (from monte carlo simulation error) for the p-values returned by a test
get.p.value <- function( tst ) {
    cnts = (tst$ci.p - tst$gamma) * tst$B
    bts = sapply( cnts, function( cnt ) {
        bt = binom.test( cnt, tst$B )
        confint( bt )[2:3]
    } )
    stopifnot( tst$p.value == max( tst$ci.p ) )
    c( p.value=max( tst$ci.p ), min.p= min( bts[1,] ), max.p=max( bts[2,] ), plug=tst$p.value.plug )
}





# Plot the curve of different p-values for different values of tau
plot.FRTCI.curve <- function( tst, true.tau=NULL, 
                              xlab=expression(tau), ylab="p-value", true.tau.col="red",
                              plot.envelope=TRUE, ci.line.col="blue", ... ) {
    
    cnts = (tst$ci.p - tst$gamma) * tst$B
    cnt <- function(x){
      bt <-  binom.test( cnts, tst$B )
      bt$conf.int
    }
    
    bts <-  sapply( cnts, cnt)
    
    plot( tst$te.vec, tst$ci.p, ylim=c(0,0.6), type="l", xlab=xlab, ylab=ylab, ...  )
    abline( v=tst$te.hat, col= ci.line.col )
    abline( v=c(tst$te.hat-tst$te.MOE,tst$te.hat+tst$te.MOE), col= ci.line.col, lty=2 )
    
    if ( plot.envelope ) {
        lines( tst$te.vec, bts[1,], lty=3, col="grey" )
        lines( tst$te.vec, bts[2,], lty=3, col="grey" )
        lines( lowess(tst$te.vec, tst$ci.p, f=1/6),lty=3,lwd=2,col="darkgray" )
        
    }
    
    rug( tst$te.vec, ticksize=0.01 )
    
    if (!is.null(true.tau) ) {
        abline( v=true.tau, lwd=2, col= true.tau.col )
        
    }
}



##
## some testing code for the functions above.
##
if  (FALSE ) {
    
    
    Z = rep( c(0,1), 100 )
    tau=4
    Y = ifelse( Z, rnorm( 100, tau), rnorm( 100, 0 ) )	
    tst = FRTCI( Y, Z, B=50)
    
    plot.FRTCI.curve( tst, tau=tau )
    tst
    
    FRTplug( Y, Z, B=200 )
    
}

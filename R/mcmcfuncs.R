library(coda)

logit <- function(p) log(p/(1-p))
expit <- function(x) 1/(1+exp(-x))

all.permutations <- function(n) {
  # all permutations
  # https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

rdirichlet <- function(n, alphavec) {
  # returns a matrix of n draws from a Dirichlet(alphavec)
  # distribution - one per row
  m <- length(alphavec)
  retval <- array(rgamma(n*m,rep(alphavec,n),1),dim=c(m,n))
  retval <- t(apply(retval, 2, function(x) x/sum(x)))
  return(retval)
}

ddirichlet <- function(ymat, alphavec, log=FALSE) {
  m <- length(alphavec)
  if(is.null(dim(ymat))) dim(ymat) <- c(1,length(ymat))
  if(m!=ncol(ymat)) stop("Dimensions of ymat must be n x m, where m is the length of alphvec")
  n <- nrow(ymat)
  const <- lgamma(sum(alphavec))-sum(lgamma(alphavec))
  retval <- const + apply( (alphavec-1)*log(t(ymat)), 2, sum)
  if(!log) retval <- exp(retval)
  return(retval)
}

rtrgeom <- function(n,eta,kmin=1,kmax=Inf) {
  u <- runif(n)
  x <- kmin-1 + log( 1-u*(1-(1-eta)^(kmax-kmin+1)) ) / log(1-eta)
  x <- pmax(kmin,ceiling(x))
  return(x)
}

dtrgeom <- function(x,eta,kmin=1,kmax=Inf,log=FALSE) {
  if(any(x<kmin | x>kmax)) stop("x values outside [kmin,kmax]")
  retval <- (x-1)*log(1-eta) + log(eta/( (1-eta)^(kmin-1)-(1-eta)^(kmax)  ))
  if(!log) retval <- exp(retval)
  return(retval)
}

#---------------------------------------------------------------------------------

make.state <- function(datlist, cpar, parlist=NULL, ...) {
  if(is.null(parlist)) parlist <- make.parlist(cpar, ...)
  state <- list(parlist=parlist,
                llike=llikef(parlist, datlist, cpar),
                lprior=lpriorf(parlist, cpar),
                accepted=TRUE)
  return(state)
}

burn.chain <- function(state, nburn=1, datlist, cpar,
                       verbose=FALSE, force=FALSE) {
  if(nburn>0) {
    for(i in 1:nburn) state <- update.state(state, datlist, cpar,
                                            verbose=verbose, force=force)
  }
  return(state)
}

run.chain <- function(state, nstore=100, nburn=0, nthin=1, datlist, cpar,
                      verbose=FALSE, force=FALSE, latent=FALSE) {
  library(coda)
  state <- burn.chain(state, nburn, datlist, cpar, verbose=verbose, force=force)
  svec <- state.as.vector(state,cpar,latent)
  smat <- array(NA,dim=c(nstore,length(svec)))
  colnames(smat) <- names(svec)
  if(verbose) {
    cat(paste0("-----\nStep ",1,"\n-----"))
    print(svec)
  }
  smat[1,] <- svec
  for(i in 2:nstore) {
    if(verbose) cat(paste0("-----\nStep ",i,"\n-----"))
    state <- burn.chain(state, nthin, datlist, cpar, verbose=verbose, force=force)
    smat[i,] <- state.as.vector(state,cpar,latent)
  }
  return(as.mcmc(smat))
}

run.chain.list <- function(nchain=1, nstore=100, nburn=0, nthin=1,
                           datlist, cpar, verbose=FALSE, force=FALSE,
                           latent=FALSE) {
  slist <- lapply(1:nchain,
                  function(i) {
                    state <- make.state(datlist, cpar)
                    smat <- run.chain(state, nstore=nstore,
                                      nburn=nburn, nthin=nthin,
                                      datlist=datlist, cpar=cpar,
                                      verbose=verbose, force=force,
                                      latent=latent)
                    return(smat)
                  })
  return(mcmc.list(slist))
}


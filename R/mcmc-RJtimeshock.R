#---------------------------------------------------------------------------------
# RJ time series with shocks

make.parlist <- function(cpar, R=1, tau.e=1,
                         ivec=NULL, svec=NULL) {
  if(is.null(ivec)) kvec <- sample(cpar$n, R, replace=TRUE)
  if(is.null(svec)) svec <- s0*rbeta(R, cpar$as, cpar$bs)

  latentnames <- NULL
  parlist <- list(R=R, tau.e=tau.e,
                  kvec=kvec, svec=svec)
  parnames <- names(parlist)
  parnames <- parnames[!(parnames%in%latentnames)]
  parlist$parnames <- parnames

  return(parlist)
}

sim.datlist <- function(n, parlist, cpar) {
  if(n!=cpar$n) stop("n must equal cpar$n")
  #kvec <- sample(cpar$n, parlist$R, replace=FALSE)
  #svec <- cpar$s0*rbeta(parlist$R, cpar$as, cpar$bs)
  x <- rep(0,n)
  x[parlist$kvec] <- parlist$svec

  y <- rnorm(n, x, 1/sqrt(parlist$tau.e))
  datlist <- list(y=y, n=n)
  return(datlist)
}

llikef <- function(parlist, datlist, cpar) {
  xvec <- rep(0,cpar$n)
  xvec[parlist$kvec] <- parlist$svec
  retval <- sum(dnorm(datlist$y, xvec, 1/sqrt(parlist$tau.e), log=TRUE))
  return(retval)
}
lpriorf <- function(parlist, cpar) {
  retval <- dgamma(parlist$tau.e, cpar$alpha, cpar$beta, log=TRUE)
  retval <- retval + sum(dbeta(parlist$svec/cpar$s0, cpar$as, cpar$bs, log=TRUE))
  retval <- retval - lchoose(cpar$n, parlist$R)
  retval <- retval + dtrgeom(parlist$R, cpar$eta, 1, cpar$Rmax)
  return(retval)
}

plot.state <- function(state, cpar, main="", ...) {
  plot(NA,NA, xlim=c(1,cpar$n), ylim=c(0,cpar$s0),
       xlab="t", ylab="s", ...)
  segments(state$parlist$kvec, rep(0,length(state$parlist$kvec)),
           state$parlist$kvec, state$parlist$svec, ...)
  points(state$parlist$kvec, state$parlist$svec, pch=16)
  title(main=main)
  invisible()
}

plot.statemat <- function(smat, cpar, thin=1) {
  ns <- nrow(smat)
  np <- floor(ns/thin)
  for(i in thin*(1:np)) {
    state <- vector.as.state(smat[i,],cpar)
    plot.state(state, cpar, main=i)
  }
  invisible()
}

set.move.prob <- function(parlist, cpar) {
  mprob <- cpar$rjmoves
  mprob["birth"] <- mprob["birth"]*(parlist$R<cpar$Rmax)
  mprob["death"] <- mprob["death"]*(parlist$R>1)
  mprob <- mprob/sum(mprob)
  return(mprob)
}

update.state <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  if(verbose) {
    cat(paste0("START: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
  }
  state <- update.state.rj(state, datlist, cpar, verbose, force)
  state <- update.state.tau.e(state, datlist, cpar, verbose, force)
  state <- update.state.kvec(state, datlist, cpar, verbose, force)
  state <- update.state.svec(state, datlist, cpar, verbose, force)

  return(state)
}

update.state.tau.e <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update tau.e: Gibbs proposal
  muvec <- rep(0,cpar$n)
  muvec[state$parlist$kvec] <- state$parlist$svec
  state$parlist$tau.e <- rgamma(1, cpar$alpha+cpar$n/2,
                                   cpar$beta + 0.5*sum( (datlist$y-muvec)^2 ))
  state$llike <- llikef(state$parlist, datlist, cpar)
  state$lprior <- lpriorf(state$parlist, cpar)
  state$accepted <- TRUE

  return(state)
}

update.state.kvec <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update kvec
  for(k in 1:state$parlist$R) {
    ostate <- state
    i1 <- state$parlist$kvec[k]
    i2 <- sample((1:cpar$n)[-state$parlist$kvec],1)
    s <- state$parlist$svec[k]
    state$parlist$kvec[k] <- i2
    log.r <- ( -0.5*state$parlist$tau.e*(
      (datlist$y[i2]-s)^2 + (datlist$y[i1])^2
     -(datlist$y[i1]-s)^2 - (datlist$y[i2])^2
    ))
    if(verbose) {
      cat(paste0("svec[",k,"]: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
      state$llike <- llikef(state$parlist, datlist, cpar)
      state$lprior <- lpriorf(state$parlist, cpar)
      state$accepted <- TRUE
    } else {
      # Reject
      state <- ostate
      state$accepted <- FALSE
    }
    if(verbose) cat(paste0("Proposal was ","not "[!state$accepted],"accepted.\n"))
  }

  return(state)
}

update.state.svec <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update svec
  for(k in 1:state$parlist$R) {
    ostate <- state
    i <- state$parlist$kvec[k]
    s1 <- state$parlist$svec[k]
    u1 <- logit(state$parlist$svec[k]/cpar$s0)
    u2 <- rnorm(1, u1, cpar$sigma.qs)
    s2 <- cpar$s0*expit(u2)
    state$parlist$svec[k] <- s2
    log.r <- ( -(0.5*state$parlist$tau.e)*( (datlist$y[i]-s2)^2 - (datlist$y[i]-s1)^2 )
               + cpar$as*log(s2/s1) + cpar$bs*log((cpar$s0-s2)/(cpar$s0-s1))
                 )
    if(verbose) {
      cat(paste0("s[",k,"]: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
      state$llike <- llikef(state$parlist, datlist, cpar)
      state$lprior <- lpriorf(state$parlist, cpar)
      state$accepted <- TRUE
    } else {
      # Reject
      state <- ostate
      state$accepted <- FALSE
    }
    if(verbose) cat(paste0("Proposal was ","not "[!state$accepted],"accepted.\n"))
  }

  return(state)
}

update.state.rj <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # NOT COMPLETE
  #if(verbose) browser() ##!!==

  # birth/death/hold
  omprob <- set.move.prob(state$parlist, cpar)
  move <- sample(3, 1, prob=omprob)
  if(move==1) {
    # Birth
    ostate <- state

    R <- state$parlist$R
    s <- s0*rbeta(1, cpar$as, cpar$bs)
    i <- sample((1:cpar$n)[-state$parlist$kvec],1)
    state$parlist$R <- R+1
    state$parlist$kvec <- c(state$parlist$kvec, i)
    state$parlist$svec <- c(state$parlist$svec, s)

    mprob <- set.move.prob(state$parlist, cpar)
    log.r <- ( -0.5*state$parlist$tau.e*( (datlist$y[i]-s)^2 - (datlist$y[i])^2 ) )
    rp <- (mprob[2]/omprob[1]) * (cpar$n-R)/(R+1) * (1-cpar$eta)

    log.r <- as.vector(log.r + log(rp))
    if(verbose) {
      cat(paste0("birth: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
      state$llike <- llikef(state$parlist, datlist, cpar)
      state$lprior <- lpriorf(state$parlist, cpar)
      state$accepted <- TRUE
    } else {
      # Reject
      state <- ostate
      state$accepted <- FALSE
    }
    if(verbose) cat(paste0("Proposal was ","not "[!state$accepted],"accepted.\n"))
  } else if(move==2) {
    # Death
    ostate <- state

    R <- state$parlist$R
    j <- sample(R,1)
    i <- state$parlist$kvec[j]
    s <- state$parlist$svec[j]
    state$parlist$R <- R-1
    state$parlist$kvec <- state$parlist$kvec[-j]
    state$parlist$svec <- state$parlist$svec[-j]

    mprob <- set.move.prob(state$parlist, cpar)
    log.r <- ( -0.5*state$parlist$tau.e*( (datlist$y[i])^2 - (datlist$y[i]-s)^2  ) )
    rp <- (mprob[1]/omprob[2]) * R/(cpar$n-R+1) / (1-cpar$eta)
    log.r <- as.vector(log.r + log(rp))
    if(verbose) {
      cat(paste0("death: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
      state$llike <- llikef(state$parlist, datlist, cpar)
      state$lprior <- lpriorf(state$parlist, cpar)
      state$accepted <- TRUE
    } else {
      # Reject
      state <- ostate
      state$accepted <- FALSE
    }
    if(verbose) cat(paste0("Proposal was ","not "[!state$accepted],"accepted.\n"))
  } else {
    if(verbose) {
      log.r <- 0
      cat(paste0("hold: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
  }

  return(state)
}


state.as.vector <- function(state, cpar, latent=FALSE) {
  if(!latent) state$parlist <- state$parlist[names(state$parlist)%in%c(state$parlist$parnames)]
  state$parlist$parnames <- NULL
  state$parlist$kvec <- c(state$parlist$kvec,rep(NA,cpar$Rmax-state$parlist$R))
  state$parlist$svec <- c(state$parlist$svec,rep(NA,cpar$Rmax-state$parlist$R))
  unlist(state)
}

vector.as.state <- function(svec, cpar) {
  parnames <- names(svec)
  parnames <- grep("^parlist.",parnames,value=TRUE)
  parvec <- svec[parnames]
  names(parvec) <- gsub("^parlist\\.","",names(parvec))
  R <- as.vector(parvec["R"])
  state <- list(parlist=list(R=R,
                             tau.e=as.vector(parvec["tau.e"]),
                             kvec=as.vector(parvec[grep("kvec",names(parvec))][1:R]),
                             svec=as.vector(parvec[grep("svec",names(parvec))][1:R]),
                             parnames=c("R","tau.e","kvec","svec")),
                llike=as.vector(svec["llike"]),
                lprior=as.vector(svec["lprior"]),
                accepted=as.logical(svec["accepted"]))
  return(state)
}

relabel.state <- function(state, neworder=NULL) {
  R <- state$parlist$R
  if(is.null(neworder)) neworder <- sample(R,R)
  state$parlist$kvec <- state$parlist$kvec[neworder]
  state$parlist$svec <- state$parlist$svec[neworder]
  return(state)
}

relabel.svec <- function(svec, neworder=NULL) {
  R <- svec["parlist.R"]
  if(is.null(neworder)) neworder <- sample(R,R)
  for(vecname in c("kvec","svec")) {
    idx <- grep(paste0("^parlist.",vecname),names(svec))[1:R]
    svec[idx] <- svec[idx][neworder]
  }
  return(svec)
}

relabel.smat <- function(smat, nmat=NULL) {
  if(is.null(nmat)) return(smat)
  smat <- t(sapply(1:nrow(smat),
                   function(im) {
                     relabel.svec(smat[im,], neworder=nmat[im,])
                   }))
  return(smat)
}


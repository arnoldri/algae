#---------------------------------------------------------------------------------
# RJ algae

make.parlist <- function(cpar,
                         alpha=1, beta=1, tau.e=1, ne=1,
                         ivec=NULL, tvec=NULL, wvec=NULL, svec=NULL) {
  if(is.null(ivec)) ivec <- sample(cpar$nm, ne, prob=cpar$pivec, replace=TRUE)
  if(is.null(tvec)) tvec <- sample(cpar$nt, ne, replace=TRUE)
  if(is.null(wvec)) wvec <- rnbinom(ne, cpar$rw, cpar$pw)
  if(is.null(svec)) svec <- rgamma(ne, cpar$as, cpar$bs)

  latentnames <- NULL
  parlist <- list(alpha=alpha, beta=beta, tau.e=tau.e, ne=ne,
                  ivec=ivec, tvec=tvec, wvec=wvec, svec=svec)
  parnames <- names(parlist)
  parnames <- parnames[!(parnames%in%latentnames)]
  parlist$parnames <- parnames

  return(parlist)
}

make.mumat <- function(parlist, cpar, .useC=TRUE) {
  if(.useC) {
    return(rcpp_make_mumat(parlist,cpar))
  } else {
    return(make.mats(parlist,cpar)$mumat)
  }
}

make.mats <- function(parlist, cpar) {
  # active matrix
  amat <- array(0,dim=c(cpar$nt+1,cpar$nm))
  if(parlist$ne>0) {
    for(k in 1:parlist$ne) {
      tsub <- parlist$tvec[k]+(0:parlist$wvec[k])
      tsub <- tsub[tsub<=cpar$nt]
      amat[1+tsub,parlist$ivec[k]] <- 1
    }
  }
  ###amat <- array(pmin(amat,1),dim=dim(amat))
  # new biomass
  bmat <- array(0,dim=c(cpar$nt+1,cpar$nm))
  if(parlist$ne>0) {
    for(k in 1:parlist$ne) {
      bmat[1+parlist$tvec[k],parlist$ivec[k]] <- parlist$svec[k]
    }
  }
  # total biomass
  mmat <- array(0,dim=c(cpar$nt+1,cpar$nm))
  if(parlist$ne>0) {
    for(t in 1+(1:cpar$nt)) {
      mmat[t,] <- ( parlist$alpha*amat[t,]*mmat[t-1,]
                      + parlist$beta*cpar$pmat%*%mmat[t-1,]
                      + bmat[t,] )
    }
  }
  mumat <- cpar$delta + mmat[-1,]/rep(cpar$volvec,each=cpar$nt)
  return(list(amat=amat, bmat=bmat, mmat=mmat, mumat=mumat))
}

sim.datlist <- function(parlist, cpar) {
  mats <- make.mats(parlist, cpar)
  ymat <- exp(rnorm(cpar$nt*cpar$nm,
                    log(mats$mumat),
                    1/sqrt(parlist$tau.e)))
  dim(ymat) <- c(cpar$nt,cpar$nm)
  datlist <- list(ymat=ymat)
  return(datlist)
}

llikef <- function(parlist, datlist, cpar) {
  mats <- make.mats(parlist, cpar)
  retval <- sum(dnorm(log(datlist$ymat),
                      log(mats$mumat), 1/sqrt(parlist$tau.e), log=TRUE),
                na.rm=TRUE)
  return(retval)
}
lpriorf <- function(parlist, cpar) {
  retval <- dtrgeom(parlist$ne, cpar$eta, 0, cpar$nemax)
  retval <- retval + dgamma(parlist$alpha, cpar$aa, cpar$ba)
  retval <- retval + dbeta(parlist$beta, cpar$ab, cpar$bb)
  retval <- retval + dgamma(parlist$tau.e, cpar$at, cpar$bt, log=TRUE)
  retval <- retval + sum(log(cpar$pivec[parlist$ivec]))
  retval <- retval - parlist$ne*log(cpar$nt)
  retval <- retval + sum(dnbinom(parlist$wvec, cpar$rw, cpar$pw, log=TRUE))
  retval <- retval + sum(dgamma(parlist$svec, cpar$as, cpar$bs, log=TRUE))
  return(retval)
}

plot.state <- function(state, datlist, cpar, main="", ...) {
  plot(NA,NA, xlim=c(cpar$tmin,cpar$tmax), ylim=c(0,max(datlist$ymat,na.rm=TRUE)),
       xlab="Time",ylab="Intensity")
  invisible(apply(datlist$ymat,2,lines,x=cpar$tmin:cpar$tmax))
  points(state$parlist$tvec, state$parlist$svec, col="red", pch=16)
  segments(state$parlist$tvec, state$parlist$svec,
           state$parlist$tvec+state$parlist$wvec, state$parlist$svec,
           lwd=2, col="red")
  if(state$parlist$ne>0) {
    text(state$parlist$tvec, state$parlist$svec, labels=state$parlist$ivec, pos=3, col="red")
  }
  mlab1 <- sprintf("(N,alpha,beta,tau,sigma)=(%d, %.3f, %.3f, %.1f, %.1f)",
                  state$parlist$ne, state$parlist$alpha, state$parlist$beta,
                  state$parlist$tau.e, 1/sqrt(state$parlist$tau.e))
  mlab2 <- sprintf("(LL,LP)=(%.2f, %.2f)", state$llike, state$lprior)
  mtext(mlab1, side=3, line=0, cex=0.7, adj=0)
  mtext(mlab2, side=3, line=0, cex=0.7, adj=1)
  title(main=main)
  invisible()
}

plot.statemat <- function(smat, datlist, cpar, thin=1) {
  ns <- nrow(smat)
  np <- floor(ns/thin)
  for(i in thin*(1:np)) {
    state <- vector.as.state(smat[i,],cpar)
    plot.state(state, datlist, cpar, main=i)
  }
  invisible()
}

set.move.prob <- function(parlist, cpar) {
  mprob <- cpar$rjmoves
  mprob["birth"] <- mprob["birth"]*(parlist$ne<cpar$nemax)
  mprob["death"] <- mprob["death"]*(parlist$ne>cpar$nemin)
  mprob <- mprob/sum(mprob)
  return(mprob)
}

##!!==

update.state <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  if(verbose) {
    cat(paste0("START: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
  }
  state <- update.state.rj(state, datlist, cpar, verbose, force)
  state <- update.state.alpha(state, datlist, cpar, verbose, force)
  state <- update.state.beta(state, datlist, cpar, verbose, force)
  state <- update.state.tau.e(state, datlist, cpar, verbose, force)
  state <- update.state.ivec.v1(state, datlist, cpar, verbose, force)
  state <- update.state.ivec.v2(state, datlist, cpar, verbose, force)
  state <- update.state.tvec.v1(state, datlist, cpar, verbose, force)
  state <- update.state.tvec.v2(state, datlist, cpar, verbose, force)
  state <- update.state.wvec(state, datlist, cpar, verbose, force)
  state <- update.state.svec(state, datlist, cpar, verbose, force)

  return(state)
}

update.state.tau.e <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update tau.e: Gibbs proposal
  mats <- make.mats(state$parlist, cpar)
  nobs <- sum(cpar$phimat[!is.na(cpar$phimat)])
  state$parlist$tau.e <- rgamma(1,
                                state$parlist$alpha+nobs/2,
                                state$parlist$beta +
                                  0.5*sum( ((log(datlist$ymat)-log(mats$mumat))^2)[!is.na(cpar$phimat)] ))
  state$llike <- llikef(state$parlist, datlist, cpar)
  state$lprior <- lpriorf(state$parlist, cpar)
  state$accepted <- TRUE

  return(state)
}

update.state.alpha <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update alpha
  ostate <- state
  a1 <- state$parlist$alpha
  u1 <- log(a1)
  u2 <- rnorm(1, u1, cpar$sigma.qa)
  a2 <- exp(u2)
  state$parlist$alpha <- a2
  state$llike <- llikef(state$parlist, datlist, cpar)
  state$lprior <- lpriorf(state$parlist, cpar)

  log.r <- ( state$llike - ostate$llike )
  log.r <- log.r + ( cpar$aa*(u2-u1) - cpar$ba*(a2-a1) )
  if(verbose) {
    cat(paste0("alpha: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
  }
  if(runif(1)<exp(log.r) || force) {
    # Accept
    state$accepted <- TRUE
  } else {
    # Reject
    state <- ostate
    state$accepted <- FALSE
  }
  if(verbose) cat(paste0("Proposal was ","not "[!state$accepted],"accepted.\n"))

  return(state)
}

update.state.beta <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update beta
  ostate <- state
  b1 <- state$parlist$beta
  u1 <- logit(b1)
  u2 <- rnorm(1, u1, cpar$sigma.qb)
  b2 <- expit(u2)
  state$parlist$beta <- b2
  state$llike <- llikef(state$parlist, datlist, cpar)
  state$lprior <- lpriorf(state$parlist, cpar)

  log.r <- ( state$llike - ostate$llike )
  log.r <- log.r + ( cpar$ab*log(b2/b1) + cpar$bb*log((1-b2)/(1-b1)) )
  if(verbose) {
    cat(paste0("beta: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
  }
  if(runif(1)<exp(log.r) || force) {
    # Accept
    state$accepted <- TRUE
  } else {
    # Reject
    state <- ostate
    state$accepted <- FALSE
  }
  if(verbose) cat(paste0("Proposal was ","not "[!state$accepted],"accepted.\n"))

  return(state)
}

update.state.ivec.v1 <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update ivec (version 1: shift to a random neighbour)

  # no update possible if there is only one location, or if there are no events
  if(cpar$nm==1) {
    state$accepted <- FALSE
    return(state)
  } else if(state$parlist$ne==0) {
    state$accepted <- TRUE
    return(state)
  }

  for(k in 1:state$parlist$ne) {
    ostate <- state
    i1 <- state$parlist$ivec[k]
    if(cpar$numnbr[i1]==0) {
      # no neighbours - no move possible
      if(verbose) {
        cat(paste0("I[",k,"](V1): No neighbours to move to")) ##!!==
      }
      state$accepted <- FALSE
      if(verbose) cat(paste0("Proposal was ","not "[!state$accepted],"accepted.\n"))
    } else {
      i2 <- sample(cpar$nlist[[k]],1)
      state$parlist$ivec[k] <- i2
      state$llike <- llikef(state$parlist, datlist, cpar)
      state$lprior <- lpriorf(state$parlist, cpar)
      log.r <- (state$llike + state$lprior - ostate$llike - ostate$lprior)
      log.r <- log.r + log(cpar$numnbr[i1]/cpar$numnbr[i2])
      if(verbose) {
        cat(paste0("I[",k,"](V1): ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
      }
      if(runif(1)<exp(log.r) || force) {
        # Accept
        state$accepted <- TRUE
      } else {
        # Reject
        state <- ostate
        state$accepted <- FALSE
      }
      if(verbose) cat(paste0("Proposal was ","not "[!state$accepted],"accepted.\n"))
    }
  }

  return(state)
}

update.state.ivec.v2 <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update ivec (version 2: shift to a random location, excluding the current one)

  # no update possible if there is only one location, or no events
  if(cpar$nm==1) {
    state$accepted <- FALSE
    return(state)
  } else if(state$parlist$ne==0) {
    state$accepted <- TRUE
    return(state)
  }

  for(k in 1:state$parlist$ne) {
    ostate <- state
    i1 <- state$parlist$ivec[k]
    i2 <- sample((1:cpar$nm)[-t1],1)
    state$parlist$ivec[k] <- i2
    state$llike <- llikef(state$parlist, datlist, cpar)
    state$lprior <- lpriorf(state$parlist, cpar)
    log.r <- (state$llike + state$lprior - ostate$llike - ostate$lprior)
    if(verbose) {
      cat(paste0("I[",k,"](V2): ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
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


update.state.tvec.v1 <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update tvec (version 1: shift by 1 up or down)

  # no update possible if there is only one time point or no events
  if(cpar$nt==1) {
    state$accepted <- FALSE
    return(state)
  } else if(state$parlist$ne==0) {
    state$accepted <- TRUE
    return(state)
  }

  for(k in 1:state$parlist$ne) {
    ostate <- state
    t1 <- state$parlist$tvec[k]
    if(t1==1) { # forced upward
      t2 <- t1+1
      q12 <- 1
    } else if(t1==cpar$nt) { # forced downward
      t2 <- t1-1
      q12 <- 1
    } else { # up or down with equal probability
      t2 <- sample(t1+c(-1,1),1)
      q12 <- 0.5
    }
    # probability of the reverse move
    if(t2==1 || t2==cpar$nt) {
      q21 <- 1
    } else {
      q21 <- 0.5
    }
    state$parlist$tvec[k] <- t2
    state$llike <- llikef(state$parlist, datlist, cpar)
    state$lprior <- lpriorf(state$parlist, cpar)
    log.r <- (state$llike + state$lprior - ostate$llike - ostate$lprior
              + log(q21/q12))
    if(verbose) {
      cat(paste0("T[",k,"](V1): ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
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

update.state.tvec.v2 <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update tvec (version 2: shift to a random time)

  # no update possible if there is only one time point or if there are no events
  if(cpar$nt==1) {
    state$accepted <- FALSE
    return(state)
  } else if(state$parlist$ne==0) {
    state$accepted <- TRUE
    return(state)
  }

  for(k in 1:state$parlist$ne) {
    ostate <- state
    t1 <- state$parlist$tvec[k]
    t2 <- sample((1:cpar$nt)[-t1],1)
    state$parlist$tvec[k] <- t2
    state$llike <- llikef(state$parlist, datlist, cpar)
    state$lprior <- lpriorf(state$parlist, cpar)
    log.r <- (state$llike + state$lprior - ostate$llike - ostate$lprior)
    if(verbose) {
      cat(paste0("T[",k,"](V2): ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
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

update.state.wvec <- function(state, datlist, cpar, verbose=FALSE, force=FALSE) {

  # update wvec (shift by 1 up or down)
  # no update possible if there are no events
  if(state$parlist$ne==0) {
    state$accepted <- TRUE
    return(state)
  }

  for(k in 1:state$parlist$ne) {
    ostate <- state
    w1 <- state$parlist$wvec[k]
    if(w1==0) { # forced upward
      w2 <- w1+1
      q12 <- 1
    } else { # up or down with equal probability
      w2 <- sample(w1+c(-1,1),1)
      q12 <- 0.5
    }
    # probability of the reverse move
    if(w2==0) {
      q21 <- 1
    } else {
      q21 <- 0.5
    }
    state$parlist$wvec[k] <- w2
    state$llike <- llikef(state$parlist, datlist, cpar)
    state$lprior <- lpriorf(state$parlist, cpar)
    log.r <- (state$llike + state$lprior - ostate$llike - ostate$lprior
              + log(q21/q12))
    if(verbose) {
      cat(paste0("W[",k,"]: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
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
  # no update possible if there are no events
  if(state$parlist$ne==0) {
    state$accepted <- TRUE
    return(state)
  }

  for(k in 1:state$parlist$ne) {
    ostate <- state
    s1 <- state$parlist$svec[k]
    u1 <- log(state$parlist$svec[k])
    u2 <- rnorm(1, u1, cpar$sigma.qs)
    s2 <- exp(u2)
    state$parlist$svec[k] <- s2
    state$llike <- llikef(state$parlist, datlist, cpar)
    state$lprior <- lpriorf(state$parlist, cpar)

    log.r <- ( state$llike - ostate$llike )
    log.r <- log.r + ( cpar$as*(u2-u1) - cpar$bs*(s2-s1) )

    if(verbose) {
      cat(paste0("s[",k,"]: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
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

  #if(verbose) browser() ##!!==

  # birth/death/hold
  omprob <- set.move.prob(state$parlist, cpar)
  move <- sample(3, 1, prob=omprob)
  if(move==1) {
    # Birth
    ostate <- state

    ne <- state$parlist$ne
    i <- sample(cpar$nm, 1, prob=cpar$pivec)
    t <- sample(cpar$nt, 1)
    w <- rnbinom(1, cpar$rw, cpar$pw)
    s <- rgamma(1, cpar$as, cpar$bs)

    state$parlist$ne <- ne+1
    state$parlist$ivec <- c(state$parlist$ivec, i)
    state$parlist$tvec <- c(state$parlist$tvec, t)
    state$parlist$wvec <- c(state$parlist$wvec, w)
    state$parlist$svec <- c(state$parlist$svec, s)

    state$llike <- llikef(state$parlist, datlist, cpar)
    state$lprior <- lpriorf(state$parlist, cpar)

    mprob <- set.move.prob(state$parlist, cpar)

    log.r <- ( state$llike - ostate$llike )
    rp <- (mprob[2]/omprob[1])*(1-cpar$eta)
    log.r <- as.vector(log.r + log(rp))

    if(verbose) {
      cat(paste0("birth: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
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

    ne <- state$parlist$ne
    k <- sample(ne,1)
    i <- state$parlist$ivec[k]
    t <- state$parlist$tvec[k]
    w <- state$parlist$wvec[k]
    s <- state$parlist$svec[k]
    state$parlist$ne <- ne-1
    state$parlist$ivec <- state$parlist$ivec[-k]
    state$parlist$tvec <- state$parlist$tvec[-k]
    state$parlist$wvec <- state$parlist$wvec[-k]
    state$parlist$svec <- state$parlist$svec[-k]

    state$llike <- llikef(state$parlist, datlist, cpar)
    state$lprior <- lpriorf(state$parlist, cpar)

    mprob <- set.move.prob(state$parlist, cpar)

    log.r <- ( state$llike - ostate$llike )
    rp <- (mprob[1]/omprob[2])/(1-cpar$eta)
    log.r <- as.vector(log.r + log(rp))

    if(verbose) {
      cat(paste0("death: ")); print(c(state.as.vector(state,cpar),log.r)) ##!!==
    }
    if(runif(1)<exp(log.r) || force) {
      # Accept
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
  state$parlist$ivec <- c(state$parlist$ivec,rep(NA,cpar$nemax-state$parlist$ne))
  state$parlist$tvec <- c(state$parlist$tvec,rep(NA,cpar$nemax-state$parlist$ne))
  state$parlist$wvec <- c(state$parlist$wvec,rep(NA,cpar$nemax-state$parlist$ne))
  state$parlist$svec <- c(state$parlist$svec,rep(NA,cpar$nemax-state$parlist$ne))
  unlist(state)
}

vector.as.state <- function(svec, cpar) {
  parnames <- names(svec)
  parnames <- grep("^parlist.",parnames,value=TRUE)
  parvec <- svec[parnames]
  names(parvec) <- gsub("^parlist\\.","",names(parvec))
  ne <- as.vector(parvec["ne"])
  state <- list(parlist=list(alpha=as.vector(parvec["alpha"]),
                             beta=as.vector(parvec["beta"]),
                             tau.e=as.vector(parvec["tau.e"]),
                             ne=ne,
                             ivec=as.vector(parvec[grep("ivec",names(parvec))][1:ne]),
                             tvec=as.vector(parvec[grep("tvec",names(parvec))][1:ne]),
                             wvec=as.vector(parvec[grep("wvec",names(parvec))][1:ne]),
                             svec=as.vector(parvec[grep("svec",names(parvec))][1:ne]),
                             parnames=c("alpha","beta","tau.e","ne",
                                        "ivec","tvec","wvec","svec")),
                llike=as.vector(svec["llike"]),
                lprior=as.vector(svec["lprior"]),
                accepted=as.logical(svec["accepted"]))
  return(state)
}

relabel.state <- function(state, neworder=NULL) {
  ne <- state$parlist$ne
  if(is.null(neworder)) neworder <- sample(ne,ne)
  state$parlist$ivec <- state$parlist$ivec[neworder]
  state$parlist$tvec <- state$parlist$tvec[neworder]
  state$parlist$wvec <- state$parlist$wvec[neworder]
  state$parlist$svec <- state$parlist$svec[neworder]
  return(state)
}

relabel.svec <- function(svec, neworder=NULL) {
  ne <- svec["parlist.ne"]
  if(is.null(neworder)) neworder <- sample(ne,ne)
  for(vecname in c("ivec","tvec","wvec","svec")) {
    idx <- grep(paste0("^parlist.",vecname),names(svec))[1:ne]
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


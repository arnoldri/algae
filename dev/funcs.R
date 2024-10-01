#########################################################################
# eDNA backtracking functions
dispfunc <- function(x,y,t, x0,y0,t0, vx,vy, ax,ay, teps) {
  # this function allows x,y to be vectors of equal length, and t to be of the same length or to be a scalar
  # it also allows x,y to be scalars and t to be a vector
  # don't mix these!
  retval <- dnorm(x,x0+vx*(t-t0),sqrt(ax*(t-t0+teps))) * dnorm(y,y0+vy*(t-t0),sqrt(ay*(t-t0+teps)))
  return(retval)
}
denratefunc <- function(x,y,t, x0,y0,t0, vx,vy, ax,ay,teps, n0=1,th=Inf,tm=Inf) {
  # this version allows x,y to be vectors of equal length
  # t must be a scalar
  if(t<t0 || t-t0>tm) {
    retval <- rep(0,length(x))
  } else {
    retval <- n0*2^(-(t-t0)/th)*dispfunc(x,y,t, x0,y0,t0, vx,vy, ax,ay, teps)
  }
  return(retval)
}
hfunc <- function(x, y, x0, y0, vx, vy, ax, ay, teps, th=Inf, tm=Inf) {
  apply(matrix(c(x,y),ncol=2), 1, function(xy) {
    integrate(denratefunc_integrand,
              lower=0, upper=tm,
              x=xy[1], y=xy[2], t=tm,
              x0=x0, y0=y0, vx=vx, vy=vy, ax=ax, ay=ay, teps=teps, th=th, tm=tm)$value
  })
}
denratefunc_integrand <- function(t0, x,y,t, x0,y0, vx,vy, ax,ay, teps, n0=1,th=Inf,tm=Inf) {
  # this version allows t0 to be a vector, (x,y,t) must be scalars
  retval <- n0*2^{-(t-t0)/th}*dispfunc(x,y,t, x0,y0,t0, vx,vy, ax,ay, teps)
  retval <- retval*ifelse(t<t0 | t-t0>tm, 0, 1)
  return(retval)
}


#########################################################################
# Binary model

# Simulate
binsimfunc <- function(nsim, parlist, fparlist) {
  nm <- fparlist$nm
  cvec <- diag(fparlist$pmat1)

  eps <- parlist$eps # Chance of spontaneous generation at a site
  lambda <- parlist$lambda # Transport survival parameter
  beta <- parlist$beta # Additional in-situ survival parameter
  mu <- parlist$mu # Chance that an existing bloom will die
  delta <- parlist$delta # Probability of a false positive (when scaled by 1-cvec)
  p <- parlist$p # Probability of a true positive ###(when scaled by cvec)

  zmat <- array(0, dim=c(nsim,nm))
  ymat <- array(0, dim=c(nsim,nm))
  qmat <- array(0, dim=c(nsim,nm))
  umat <- array(0, dim=c(nsim,nm))
  xmat <- array(1, dim=c(nsim,nm))

  qmat[1,] <- eps*cvec
  for(t in 2:nsim) {
    qmat[t,] <- pmin(1,eps*cvec + beta*zmat[t-1,] + lambda*pmat1%*%zmat[t-1,])
    umat[t,] <- rbinom(nm, 1, mu)
    xmat[t,] <- 1-(pmat1%*%umat[t,]>0)
    zmat[t,] <- rbinom(nm, 1, xmat[t,]*qmat[t,])
    ymat[t,] <- rbinom(nm, 1, pmin(1,delta*(1-cvec) + p*zmat[t,]))
  }
  return(list(qmat=qmat, umat=umat, xmat=xmat, zmat=zmat, ymat=ymat))
}

# Model with residency

#' Indicator of season (off = 0 or on = 1)
#' On season is December-July
#' Off season is August-November
#'
#' @param week Week number or a vector of week numbers (Week 2505 starts 2018-01-01)
#'
#' @return 1 if the week is in the December-July on season, otherwise 0
#'
#' @export
seasonfunc <- function(week) {
  as.numeric(format(week.as.date(week), "%m")) %in% c(12,1,2,3,4,5,6,7)
}

#' Make a step function from vector vvec:
#' Insert zero values when there is a transition
#' from a zero to one, or the reverse.
#'
#' @export
block01 <- function(tvec,vvec,addheadtail=TRUE) {
  aa <- cbind(tvec,vvec)
  idxp1 <- which(diff(aa[,2])==+1)
  idxm1 <- which(diff(aa[,2])==-1)
  aa <- as.list(as.data.frame(t(aa)))
  for(i in idxp1) {
    aa[[i+1]] <- c(aa[[i+1]][1],0 ,aa[[i+1]])
  }
  for(i in idxm1) {
    aa[[i]] <- c(aa[[i]], aa[[i]][1],0)
  }
  aa <- t(matrix(unlist(aa),nrow=2))
  if(addheadtail) {
    aa <- rbind(c(aa[1,1],0),
                aa,
                c(aa[nrow(aa)],0))
  }
  return(aa)
}

############################################################
#' @export
find.events <- function(xmat, quantile.threshold=0.75) {
  # find locations where there is a peak in the maximum value of the
  # rows in the matrix xmat
  nt <- nrow(xmat)
  xx <- apply(xmat,1,max)
  uq <- quantile(xx,probs=quantile.threshold)
  idx1 <- (1:nt)[xx>uq]
  dd <- diff(idx1)
  idx2a <- idx1[c(Inf,dd)>1]
  idx2b <- idx1[c(dd,Inf)>1]
  xmax <- apply(cbind(idx2a,idx2b),1, function(ii) max(xx[ii[1]:ii[2]]))
  data.frame(event=1:length(idx2a),
             istart=idx2a, iend=idx2b,
             xmax=xmax)
}
############################################################
sfreset <- function() invisible(dev.off())

plot.sobj <- function(sobj, type="polygons") {
  if(type=="polygons") {
    plot((sobj$polysf %>% st_geometry()), main=arealabel,reset=FALSE)
    text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
         lab=as.character((sobj$polysf %>% st_drop_geometry())[,sobj$clusteridname]), cex=0.5, col="black")
  } else if(type=="volume") {
    plot((sobj$polysf %>% select(volume_15)), main="Polygon Volumes to 15m")
  } else if(type=="retention") {
    plot((sobj$polysf %>% select(retention)), main="Polygon 1 week retention")
  } else if(type=="residency") {
    plot((sobj$polysf %>% select(residency)), main="Polygon 1 week residency")
  } else if(type=="communication") {
    image(1:sobj$npoly, 1:sobj$npoly, sobj$pmat, axes=FALSE, xlab="Destination", ylab="Source",
          main="Communication between polygons",asp=1)
    axis(1,at=1:sobj$npoly,lab=sobj$polyid,las=2,cex.axis=0.7)
    axis(2,at=1:sobj$npoly,lab=sobj$polyid,las=2,cex.axis=0.7)
    box()
  } else {
    invisible()
  }
  invisible()
}
do.sim <- function(sobj, tlist, parlist) {
  simlist <- list()

  # Epochs
  evec <- rep(0,tlist$nt)
  evec[1] <- rbinom(1,1,parlist$tau0*tlist$seasonvec[1])
  for(i in 2:tlist$nt) {
    evec[i] <- if(evec[i-1]==0) rbinom(1,1,parlist$tau0*tlist$seasonvec[i]) else rbinom(1,1,parlist$tau1*tlist$seasonvec[i])
  }
  simlist$evec <- evec

  # Bloom initiation probability
  pavec <- sobj$polysf$volume_15*exp(parlist$alpha0 + parlist$lambda0*log(sobj$polysf$residency))
  simlist$pavec <- pavec

  # Time dependent bloom initiation probabilities
  pimat <- outer(seasonvec*evec,pavec)
  colnames(pimat) <- names(pavec)
  # Innovation indicator
  imat <- array(rbinom(length(pimat),1,pimat),dim=dim(pimat))
  colnames(imat) <- names(pavec)
  simlist$pimat <- pimat
  simlist$imat <- imat

  # Biomass innovation
  btmat <- array(rgamma(prod(dim(pimat)),abio,bbio),dim=dim(pimat))
  btmat <- imat*btmat
  simlist$btmat <- btmat

  # Mass evolution
  # growth matrix
  qmat <- array(NA,dim=dim(pimat))
  # transported mass matrix
  dmat <- array(NA,dim=dim(pimat))
  # final mass matrix
  mmat <- array(0,dim=dim(pimat))
  # final density matrix
  cmat <- mmat
  # initiation
  mmat[1,] <- btmat[1,]
  cmat[1,] <- mmat[1,]/volvec
  qmat[1,] <- 1
  dmat[1,] <- 0

  volvec <- sobj$polysf$volume_15
  resvec <- sobj$resvec
  rhocvec <- exp(parlist$alphac + parlist$lambdac*log(sobj$resvec))
  mcvec <- volvec*rhocvec
  # post-transport variation
  etamat <- array(rnorm(prod(dim(pimat)),0,parlist$sigma.eta),dim=dim(pimat))
  for(t in 2:nt) {
    # growth rate
    qvec <- exp(-parlist$gamma2*(1-evec[t]) + evec[t]*(parlist$alpha2+parlist$lambda2*log(sobj$resvec)))
    # pre-transport growth
    gvec <- exp(evec[t]*parlist$lambda1*log(sobj$resvec))
    # transported mass
    dvec <- btmat[t,] + sobj$pmat%*%diag(gvec)%*%mmat[t-1,]
    # final grown mass with post-transport variation
    mvec <- exp(etamat[t,])*mcvec*(1-exp(-qvec*dvec/mcvec))
    # store
    mmat[t,] <- mvec
    cmat[t,] <- mmat[t,]/volvec
    dmat[t,] <- dvec
    qmat[t,] <- qvec
  }
  # observed concentration
  ymat <- array(exp(rnorm(prod(dim(cmat)),0,sigma.logy)),dim=dim(pimat))
  simlist <- c(simlist,list(mmat=mmat,cmat=cmat,dmat=dmat,qmat=qmat))

  return(simlist)
}
plot.sim <- function(simlist, sobj, tlist, parlist, type="epoch") {
  if(type=="epoch") {
    plot(week.as.date(tlist$tvec), tlist$seasonvec, type="s", col="grey",
         xlab="Date", ylab="Bloom favourability",
         main="Bloom epochs", axes=FALSE)
    bb <- block01(tlist$tvec,tlist$seasonvec)
    polygon(week.as.date(bb[,1]), bb[,2], col="grey")
    lines(week.as.date(tlist$tvec), simlist$evec, type="s")
    bb <- block01(tlist$tvec,simlist$evec)
    polygon(week.as.date(bb[,1]),bb[,2], col="orange")
    #polygon(week.as.date(c(tvec[1],tvec,tvec[length(tvec)])), c(0,evec,0), col="orange")
    axis.Date(1); axis(2,at=c(0,1)); box()
  } else if(type=="carrying.density") {
    volvec <- sobj$polysf$volume_15
    rhocvec <- exp(parlist$alphac + parlist$lambdac*log(sobj$resvec))
    plot(sobj$resvec, rhocvec, xlab="Residency", ylab=expression("Polygon carrying density: "*rho[cit]),
         main="Polygon carrying density", pch=" "); text(resvec, rhocvec, lab=names(resvec), cex=0.6)
  } else if(type=="map.carrying.density") {
    volvec <- sobj$polysf$volume_15
    rhocvec <- exp(parlist$alphac + parlist$lambdac*log(sobj$resvec))
    plot((sobj$polysf %>% mutate(rhoc=rhocvec))['rhoc'], main="Carrying density",reset=FALSE)
    text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
         lab=as.character((sobj$polysf %>% st_drop_geometry())[,sobj$clusteridname]), cex=0.6, col="white")
  } else if(type=="carrying.mass") {
    volvec <- sobj$polysf$volume_15
    rhocvec <- exp(parlist$alphac + parlist$lambdac*log(sobj$resvec))
    mcvec <- volvec*rhocvec
    plot(sobj$resvec, mcvec, xlab="Residency", ylab=expression("Polygon carrying density: "*M[cit]),
         main="Polygon carrying mass", pch=" ");
    text(sobj$resvec, mcvec, lab=names(sobj$resvec), cex=0.6)
  } else if(type=="map.carrying.mass") {
    volvec <- sobj$polysf$volume_15
    rhocvec <- exp(parlist$alphac + parlist$lambdac*log(sobj$resvec))
    mcvec <- volvec*rhocvec
    plot((sobj$polysf %>% mutate(mc=mcvec))['mc'], main="Carrying mass",reset=FALSE)
    text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
         lab=as.character((sobj$polysf %>% st_drop_geometry())[,sobj$clusteridname]), cex=0.6, col="white")
  } else if(type=="bloom.initiation.probability") {
    pavec <- sobj$polysf$volume_15*exp(parlist$alpha0 + parlist$lambda0*log(sobj$resvec))
    plot(sobj$resvec, pavec, xlab="Residency", ylab=expression("Bloom initiation probability: "*pi[it]),
         main="Bloom initiation probability", pch=" ");
    text(sobj$resvec, pavec, lab=names(sobj$resvec), cex=0.6)
  } else if(type=="map.bloom.initiation.probability") {
    pavec <- sobj$polysf$volume_15*exp(parlist$alpha0 + parlist$lambda0*log(sobj$resvec))
    plot((sobj$polysf %>% mutate(pa=pavec))['pa'], main="Bloom initiation probability", reset=FALSE)
    text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
         lab=as.character((sobj$polysf %>% st_drop_geometry())[,sobj$clusteridname]), cex=0.6, col="white")
  } else if(type=="map.number.of.initiations") {
    nivec <- apply(simlist$imat,2,sum)
    breaks <- -0.5+0:(1+max(nivec))
    plot((sobj$polysf %>% mutate(nstart=nivec))['nstart'],
         main="Number of initiations", reset=FALSE, breaks=breaks)
    text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
         lab=as.character((sobj$polysf %>% st_drop_geometry())[,sobj$clusteridname]), cex=0.6, col="white")
  } else if(type=="bloom.initiation.times") {
    image(week.as.date(tlist$tvec), 1:sobj$npoly, simlist$imat, axes=FALSE, xlab="Week", ylab="Polygon")
    title("Bloom initiation")
    axis.Date(1); axis(2,at=1:n,lab=names(simlist$pavec), las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)
  } else if(type=="innovation.biomass") {
    image(week.as.date(tlist$tvec), 1:sobj$npoly, simlist$btmat, axes=FALSE, xlab="Week", ylab="Polygon")
    title("Innovation biomass")
    axis.Date(1,format="%Y"); axis(2,at=1:n,names(simlist$pavec),las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)
  } else if(type=="biomass.over.time") {
    image(week.as.date(tlist$tvec), 1:sobj$npoly, simlist$mmat, axes=FALSE, xlab="Week", ylab="Polygon")
    title("Algal biomass")
    axis.Date(1,format="%Y"); axis(2,at=1:sobj$npoly,names(simlist$pavec),las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)
  } else if(type=="concentration.over.time") {
    image(week.as.date(tvec), 1:n, cmat, axes=FALSE, xlab="Week", ylab="Polygon")
    title("Algal concentration")
    axis.Date(1,format="%Y"); axis(2,at=1:sobj$npoly,names(simlist$pavec),las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)
  } else if(type=="observed.concentration.over.time") {
    image(week.as.date(tvec), 1:n, cmat, axes=FALSE, xlab="Week", ylab="Polygon")
    title("Observed algal concentration")
    axis.Date(1,format="%Y"); axis(2,at=1:sobj$npoly,names(simlist$pavec),las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)
  } else {
    invisible()
  }
  invisible()
}


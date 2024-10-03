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
  xx <- apply(xmat,1,max,na.rm=TRUE)
  uq <- quantile(xx,probs=quantile.threshold,na.rm=TRUE)
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
#' Reset the current device - needed after some mapping images
#'
#' @export
sfreset <- function() invisible(dev.off())

#' @export
make.datlist <- function(sobj, eventdat, monitorid) {
  datlist <- list(event=eventdat, monitorid=monitorid)

  week1 <- min(eventdat$week)
  week2 <- max(eventdat$week)
  tvec <- week1:week2
  nt <- length(tvec)
  datevec <- week.as.date(tvec)
  year1 <-as.numeric(format(datevec[1],"%Y"))
  year2 <-as.numeric(format(datevec[nt],"%Y"))

  polyid <- sort(unique(eventdat$polyid))
  cpolyid <- as.character(polyid)
  npoly <- length(polyid)
  poly.idx <- match(polyid,sobj$polysf$polyid)

  cmat <- array(NA,dim=c(nt,npoly))
  dimnames(cmat) <- list(tvec, polyid)
  cmat[cbind(match(eventdat$week, tvec),
             match(eventdat$polyid, colnames(cmat)))] <- eventdat$Value
  ievents <- find.events(cmat)
  ievents$weekstart <- tvec[ievents$istart]
  ievents$weekend <- tvec[ievents$iend]

  datlist <- list(event=eventdat, monitorid=monitorid,
                  nt=nt, tvec=tvec, datevec=datevec,
                  week1=week1, year1=year1, year2=year2,
                  npoly=npoly, polyid=polyid, cpolyid=cpolyid,
                  poly.idx=poly.idx,
                  cmat=cmat, ievents=ievents)
  return(datlist)
}

plot.sobj <- function(sobj, datlist, type="polygons",
                      main=NULL, monitored=FALSE,
                      fill.coast=TRUE, coast.col="dark green") {

  if(type=="polygons") {
    if(is.null(main)) main <- regionlabel
    if(monitored) {
      plot((sobj$polysf %>% filter(polyid%in%datlist$monitorid) %>%
              st_geometry()), main=main,reset=FALSE)
      text(sobj$polysf %>% filter(polyid%in%datlist$monitorid) %>%
             st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
           lab=unlist((sobj$polysf %>% st_drop_geometry())['polyid']), cex=0.5, col="black")
      mtext("Monitored Polygons", side=3, line=0, adj=0, cex=0.8)
    } else {
      plot((sobj$polysf %>% st_geometry()), main=regionlabel,reset=FALSE)
      text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
           lab=unlist((sobj$polysf %>% st_drop_geometry())['polyid']), cex=0.5, col="black")
    }
    if(fill.coast) plot(st_geometry(sobj$coastsf), col=coast.col, add=TRUE)

  } else if(type=="volume") {
    if(is.null(main)) main <- "Polygon Volumes to 15m"
    plot((sobj$polysf %>% select(volume_15)), main=main, reset=FALSE)
    if(fill.coast) plot(st_geometry(sobj$coastsf), col=coast.col, add=TRUE)

  } else if(type=="retention") {
    if(is.null(main)) main <- "Polygon 1 week retention"
    plot((sobj$polysf %>% select(retention)), main=main, reset=FALSE)
    if(fill.coast) plot(st_geometry(sobj$coastsf), col=coast.col, add=TRUE)

  } else if(type=="residency") {
    if(is.null(main)) main <- "Polygon 1 week residency"
    plot((sobj$polysf %>% select(residency)), main=main, reset=FALSE)
    if(fill.coast) plot(st_geometry(sobj$coastsf), col=coast.col, add=TRUE)

  } else if(type=="communication") {
    if(is.null(main)) main <- "Communication between polygons"
    image(1:sobj$npoly, 1:sobj$npoly, sobj$pmat, axes=FALSE, xlab="Destination", ylab="Source",
          main=main,asp=1)
    axis(1,at=1:sobj$npoly,lab=sobj$polyid,las=2,cex.axis=0.7)
    axis(2,at=1:sobj$npoly,lab=sobj$polyid,las=2,cex.axis=0.7)
    box()

  } else {
    warning(paste0("Plot type ",type," not recognised."))
  }
  invisible()
}

do.sim <- function(sobj, tlist, parlist, datlist) {
  simlist <- list()

  # Volume and residency
  volvec <- sobj$volvec
  resvec <- sobj$resvec

  # Epochs
  evec <- rep(0,tlist$nt)
  evec[1] <- rbinom(1,1,parlist$tau0*tlist$seasonvec[1])
  for(it in 2:tlist$nt) {
    evec[it] <- if(evec[it-1]==0) rbinom(1,1,parlist$tau0*tlist$seasonvec[it]) else rbinom(1,1,parlist$tau1*tlist$seasonvec[it])
  }
  simlist$evec <- evec

  # Bloom initiation probability
  pavec <- volvec*exp(parlist$alpha0 + parlist$lambda0*log(resvec))
  simlist$pavec <- pavec

  # Time dependent bloom initiation probabilities
  pimat <- outer(seasonvec*evec,pavec)
  dimnames(pimat) <- list(tlist$tvec, sobj$polyid)
  # Innovation indicator
  imat <- array(rbinom(length(pimat),1,pimat),dim=dim(pimat))
  dimnames(imat) <- dimnames(pimat)
  simlist$pimat <- pimat
  simlist$imat <- imat

  # Biomass innovation
  btmat <- array(rgamma(prod(dim(pimat)),abio,bbio),dim=dim(pimat))
  btmat <- imat*btmat
  dimnames(btmat) <- dimnames(pimat)
  simlist$btmat <- btmat

  # Mass evolution
  # growth matrix
  qmat <- array(NA,dim=dim(pimat),dimnames=dimnames(pimat))
  # transported mass matrix
  dmat <- array(NA,dim=dim(pimat),dimnames=dimnames(pimat))
  # final mass matrix
  mmat <- array(0,dim=dim(pimat),dimnames=dimnames(pimat))
  # final density matrix
  cmat <- mmat
  # initiation
  mmat[1,] <- btmat[1,]
  cmat[1,] <- mmat[1,]/sobj$volvec
  qmat[1,] <- 1
  dmat[1,] <- 0

  rhocvec <- exp(parlist$alphac + parlist$lambdac*log(resvec))
  mcvec <- volvec*rhocvec
  # post-transport variation
  etamat <- array(rnorm(prod(dim(pimat)),0,parlist$sigma.eta),dim=dim(pimat))
  for(it in 2:nt) {
    # growth rate
    qvec <- exp(-parlist$gamma2*(1-evec[it]) + evec[it]*(parlist$alpha2+parlist$lambda2*log(sobj$resvec)))
    # pre-transport growth
    gvec <- exp(evec[it]*parlist$lambda1*log(sobj$resvec))
    # transported mass
    dvec <- btmat[it,] + sobj$pmat%*%diag(gvec)%*%mmat[it-1,]
    # final grown mass with post-transport variation
    mvec <- exp(etamat[it,])*mcvec*(1-exp(-qvec*dvec/mcvec))
    # store
    mmat[it,] <- mvec
    cmat[it,] <- mmat[it,]/volvec
    dmat[it,] <- dvec
    qmat[it,] <- qvec
  }
  ievents <- find.events(cmat)
  ievents$weekstart <- tlist$tvec[ievents$istart]
  ievents$weekend <- tlist$tvec[ievents$iend]

  # observed concentration
  ymat <- array(exp(rnorm(prod(dim(cmat)),log(cmat),sigma.logy)),dim=dim(pimat))
  ymat[cmat==0] <- 0
  dimnames(ymat) <- dimnames(cmat)
  simlist <- c(simlist,list(mmat=mmat,cmat=cmat,dmat=dmat,qmat=qmat,
                            ymat=ymat,ievents=ievents))
  return(simlist)
}

plot.sim <- function(simlist, sobj, tlist, parlist, datlist, type="epoch",
                     main=NULL, monitored=FALSE,
                     logscale=FALSE, logoffset=1, date.axis=TRUE,
                     fill.coast=TRUE, coast.col="dark green",
                     colvec=NULL, crange=NULL, it1=NULL, it2=NULL, w1=NULL, w2=NULL,
                     make.movie=FALSE, movie.name="movie", movie.rootdir="ignore/",
                     movie.fps=20, movie.clean=TRUE) {
  if(type=="epoch") {
    if(is.null(main)) main <- "Bloom epochs"
    plot(week.as.date(tlist$tvec), tlist$seasonvec, type="s", col="grey",
         xlab="Date", ylab="Bloom favourability",
         main=main, axes=FALSE)
    bb <- block01(tlist$tvec,tlist$seasonvec)
    polygon(week.as.date(bb[,1]), bb[,2], col="grey")
    lines(week.as.date(tlist$tvec), simlist$evec, type="s")
    bb <- block01(tlist$tvec,simlist$evec)
    polygon(week.as.date(bb[,1]),bb[,2], col="orange")
    axis.Date(1); axis(2,at=c(0,1)); box()

  } else if(type=="carrying.density") {
    if(is.null(main)) main <- "Polygon carrying density"
    volvec <- sobj$polysf$volume_15
    rhocvec <- exp(parlist$alphac + parlist$lambdac*log(sobj$resvec))
    plot(sobj$resvec, rhocvec, xlab="Residency", ylab=expression("Polygon carrying density: "*rho[cit]),
         main=main, pch=" "); text(resvec, rhocvec, lab=names(resvec), cex=0.6)

  } else if(type=="map.carrying.density") {
    if(is.null(main)) main <- "Carrying density"
    volvec <- sobj$polysf$volume_15
    rhocvec <- exp(parlist$alphac + parlist$lambdac*log(sobj$resvec))
    plot((sobj$polysf %>% mutate(rhoc=rhocvec))['rhoc'], main=main,reset=FALSE)
    text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
         lab=unlist((sobj$polysf %>% st_drop_geometry())['polyid']), cex=0.6, col="white")

  } else if(type=="carrying.mass") {
    if(is.null(main)) main <- "Polygon carrying mass"
    volvec <- sobj$polysf$volume_15
    rhocvec <- exp(parlist$alphac + parlist$lambdac*log(sobj$resvec))
    mcvec <- volvec*rhocvec
    plot(sobj$resvec, mcvec, xlab="Residency", ylab=expression("Polygon carrying density: "*M[cit]),
         main=main, pch=" ");
    text(sobj$resvec, mcvec, lab=names(sobj$resvec), cex=0.6)

  } else if(type=="map.carrying.mass") {
    if(is.null(main)) main <- "Carrying mass"
    volvec <- sobj$polysf$volume_15
    rhocvec <- exp(parlist$alphac + parlist$lambdac*log(sobj$resvec))
    mcvec <- volvec*rhocvec
    plot((sobj$polysf %>% mutate(mc=mcvec))['mc'], main=main,reset=FALSE)
    text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
         lab=unlist((sobj$polysf %>% st_drop_geometry())['polyid']), cex=0.6, col="white")

  } else if(type=="bloom.initiation.probability") {
    if(is.null(main)) main <- "Bloom initiation probability"
    pavec <- sobj$polysf$volume_15*exp(parlist$alpha0 + parlist$lambda0*log(sobj$resvec))
    plot(sobj$resvec, pavec, xlab="Residency", ylab=expression("Bloom initiation probability: "*pi[it]),
         main=main, pch=" ");
    text(sobj$resvec, pavec, lab=names(sobj$resvec), cex=0.6)

  } else if(type=="map.bloom.initiation.probability") {
    if(is.null(main)) main <- "Bloom initiation probability"
    pavec <- sobj$polysf$volume_15*exp(parlist$alpha0 + parlist$lambda0*log(sobj$resvec))
    plot((sobj$polysf %>% mutate(pa=pavec))['pa'], main=main, reset=FALSE)
    text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
         lab=unlist((sobj$polysf %>% st_drop_geometry())['polyid']), cex=0.6, col="white")

  } else if(type=="map.number.of.initiations") {
    if(is.null(main)) main <- "Number of initiations"
    nivec <- apply(simlist$imat,2,sum)
    breaks <- -0.5+0:(1+max(nivec))
    plot((sobj$polysf %>% mutate(nstart=nivec))['nstart'],
         main=main, reset=FALSE, breaks=breaks)
    text(sobj$polysf %>% st_centroid(of_largest_polygon=TRUE) %>% st_coordinates(),
         lab=unlist((sobj$polysf %>% st_drop_geometry())['polyid']), cex=0.6, col="white")

  } else if(type=="bloom.initiation.times") {
    if(is.null(main)) main <- "Bloom initiation"
    image(week.as.date(tlist$tvec), 1:sobj$npoly, simlist$imat, axes=FALSE, xlab="Week", ylab="Polygon")
    title(main)
    axis.Date(1); axis(2,at=1:sobj$npoly,lab=names(simlist$pavec), las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)

  } else if(type=="innovation.biomass") {
    if(is.null(main)) main <- "Innovation biomass"
    image(week.as.date(tlist$tvec), 1:sobj$npoly, simlist$btmat, axes=FALSE, xlab="Week", ylab="Polygon")
    title(main)
    axis.Date(1,format="%Y"); axis(2,at=1:sobj$npoly,names(simlist$pavec),las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)

  } else if(type=="biomass.over.time") {
    if(is.null(main)) main <- "Algal biomass"
    if(monitored) {
      idx <- match(datlist$monitorid,sobj$polyid)
      pid <- datlist$monitorid
    } else {
      idx <- 1:sobj$npoly
      pid <- sobj$polyid
    }
    nplot <- length(idx)
    image(week.as.date(tlist$tvec), 1:nplot, simlist$mmat[,idx], axes=FALSE, xlab="Week", ylab="Polygon")
    title(main)
    axis.Date(1,format="%Y"); axis(2,at=1:nplot,as.character(pid),las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)

  } else if(type=="concentration.over.time") {
    if(monitored) {
      idx <- match(datlist$monitorid,sobj$polyid)
      pid <- datlist$monitorid
    } else {
      idx <- 1:sobj$npoly
      pid <- sobj$polyid
    }
    nplot <- length(idx)
    if(logscale) {
      if(is.null(main)) main <- "Algal concentration (log scale)"
      image(week.as.date(tlist$tvec), 1:nplot, log(logoffset+simlist$cmat[,idx]), axes=FALSE, xlab="Week", ylab="Polygon")
      title(main)
    } else {
      if(is.null(main)) main <- "Algal concentration"
      image(week.as.date(tlist$tvec), 1:nplot, simlist$cmat[,idx], axes=FALSE, xlab="Week", ylab="Polygon")
      title(main)
    }
    axis.Date(1,format="%Y"); axis(2,at=1:nplot,as.character(pid),las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)

  } else if(type=="observed.concentration.over.time") {
    if(monitored) {
      idx <- match(datlist$monitorid,sobj$polyid)
      pid <- datlist$monitorid
    } else {
      idx <- 1:sobj$npoly
      pid <- sobj$polyid
    }
    nplot <- length(idx)
    if(logscale) {
      if(is.null(main)) main <- "Observed algal concentration (log scale)"
      image(week.as.date(tlist$tvec), 1:nplot, log(logoffset+simlist$ymat[,idx]), axes=FALSE, xlab="Week", ylab="Polygon")
      title(main)
    } else {
      if(is.null(main)) main <- "Observed algal concentration"
      image(week.as.date(tlist$tvec), 1:nplot, simlist$ymat[,idx], axes=FALSE, xlab="Week", ylab="Polygon")
      title(main)
    }
    axis.Date(1,format="%Y"); axis(2,at=1:nplot,as.character(pid),las=2); box()
    abline(v=as.Date(paste(tlist$year1:tlist$year2,"01-01",sep="-")), lty=2)

  } else if(type=="map.concentration") {
    if(monitored) {
      idx <- match(datlist$monitorid,sobj$polyid)
      pid <- datlist$monitorid
    } else {
      idx <- 1:sobj$npoly
      pid <- sobj$polyid
    }
    nplot <- length(idx)
    if(is.null(it1)) it1 <- which(tlist$tvec==w1)
    if(is.null(it2)) {
      if(is.null(w2)) it2 <- it1 else it2 <- which(tlist$tvec==w2)
    }
    if(is.null(colvec)) {
       if(is.null(crange)) {
         if(logscale) {
           crange <- range(log(logoffset+simlist$cmat[,idx]),na.rm=TRUE)
           if(diff(crange)==0) crange <- crange[1]+c(0,1)
         } else {
           crange <- c(0,max(simlist$cmat[,idx],na.rm=TRUE))
           if(diff(crange)==0) crange <- c(0,1)
         }
       }
       nbreaks <- 11
       breaks <- seq(from=crange[1], to=crange[2], length=nbreaks)
       colvec <- colorRampPalette(c("light grey","red"))(nbreaks-1)
    }
    if(make.movie) {
      outdir <- paste0(movie.rootdir,"/",movie.name,"/")
      if(!dir.exists(outdir)) dir.create(outdir)
      par(opar); par(oma=c(3.1,0.5,0,0)); par(mar=1.1*c(2.5,1,4.1,4.1))
      fstem <- movie.name
      fmovie <- paste0(movie.rootdir,"/",movie.name,".gif")
    }
    for(it in it1:it2) {
      if(make.movie) {
         fname <- sprintf("%s/%s%04d.png", outdir, fstem, it)
         png(file=fname, width=480, height=480)
         par(oma=c(3.1,0.5,0,0))
         par(mar=1.1*c(2.5,1,4.1,4.1))
       }
       week <- tlist$tvec[it]
       if(is.null(main)) main <- paste0("Algal concentration in Week ",week)
       if(logscale) {
         plot((sobj$polysf %>% mutate(concentration=log(logoffset+simlist$cmat[it,])) %>% filter(polyid%in%pid))['concentration'],
            pal=colvec,breaks=breaks,
            main=main,reset=FALSE)
       } else {
         plot((sobj$polysf %>% mutate(concentration=simlist$cmat[it,]) %>% filter(polyid%in%pid))['concentration'],
              pal=colvec,breaks=breaks,
              main=main,reset=FALSE)
       }
       if(fill.coast) plot(st_geometry(sobj$coastsf), col=coast.col, add=TRUE)
       mtext(ifelse(logscale,"Log scale","Linear scale"), side=3, line=0, adj=0, cex=0.8)
       if(make.movie) {
         box()
         if(date.axis) {
           par(usr=c(as.numeric(week.as.date(tlist$tvec[c(it1,it2)])),0,1))
           axis.Date(1)
           points(week.as.date(week),0,pch=16,cex=2,xpd=TRUE)
         }
         mtext(format(week.as.date(week), "%d-%m-%Y"), side=3, line=0, adj=1, cex=0.8)
         dev.off()
       }
    }
    if(make.movie) {
      image_list <- lapply(list.files(outdir,full.names=TRUE), image_read)
      poly_animated <- image_animate(image_join(image_list), fps=movie.fps)
      image_write(poly_animated, fmovie)
      if(movie.clean) {
         unlink(list.files(outdir,full.names=TRUE))
         unlink(outdir, recursive=TRUE)
      }
    }

  } else {
    warning(paste0("Plot type ",type," not recognised."))
  }
  invisible()
}


plot.data <- function(sobj, datlist, type="epoch", main=NULL,
                      logscale=FALSE, logoffset=1, date.axis=TRUE,
                      fill.coast=TRUE, coast.col="dark green",
                      colvec=NULL, crange=NULL, it1=NULL, it2=NULL, w1=NULL, w2=NULL,
                      make.movie=FALSE, movie.name="datamovie", movie.rootdir="ignore/",
                      movie.fps=20, movie.clean=TRUE) {

  if(type=="observed.concentration.over.time") {
    if(logscale) {
      if(is.null(main)) main <- "Observed algal concentration (log scale)"
      image(week.as.date(datlist$tvec), 1:datlist$npoly, log(logoffset+datlist$cmat), axes=FALSE, xlab="Week", ylab="Polygon")
      title(main)
    } else {
      if(is.null(main)) main <- "Observed algal concentration"
      image(week.as.date(datlist$tvec), 1:datlist$npoly, datlist$cmat, axes=FALSE, xlab="Week", ylab="Polygon")
      title(main)
    }
    axis.Date(1,format="%Y"); axis(2,at=1:datlist$npoly,datlist$polyid,las=2); box()
    abline(v=as.Date(paste(datlist$year1:datlist$year2,"01-01",sep="-")), lty=2)

  } else if(type=="map.concentration") {

    if(is.null(it1)) it1 <- which(datlist$tvec==w1)
    if(is.null(it2)) {
      if(is.null(w2)) it2 <- it1 else it2 <- which(datlist$tvec==w2)
    }
    if(is.null(colvec)) {
      if(is.null(crange)) {
        if(logscale) {
          crange <- range(log(logoffset+datlist$cmat),na.rm=TRUE)
          if(diff(crange)==0) crange <- crange[1]+c(0,1)
        } else {
          crange <- c(0,max(datlist$cmat,na.rm=TRUE))
          if(diff(crange)==0) crange <- c(0,1)
        }
      }
      nbreaks <- 11
      breaks <- seq(from=crange[1], to=crange[2], length=nbreaks)
      colvec <- colorRampPalette(c("light grey","red"))(nbreaks-1)
    }
    if(make.movie) {
      outdir <- paste0(movie.rootdir,"/",movie.name,"/")
      if(!dir.exists(outdir)) dir.create(outdir)
      par(opar); par(oma=c(3.1,0.5,0,0)); par(mar=1.1*c(2.5,1,4.1,4.1))
      fstem <- movie.name
      fmovie <- paste0(movie.rootdir,"/",movie.name,".gif")
    }
    for(it in it1:it2) {
      if(make.movie) {
        fname <- sprintf("%s/%s%04d.png", outdir, fstem, it)
        png(file=fname, width=480, height=480)
        par(oma=c(3.1,0.5,0,0))
        par(mar=1.1*c(2.5,1,4.1,4.1))
      }
      week <- tlist$tvec[it]
      if(is.null(main)) main <- paste0("Observed algal concentration in Week ",week)
      if(logscale) {
        plot((sobj$polysf[datlist$poly.idx,] %>%
              mutate(concentration=log(logoffset+datlist$cmat[it,])))['concentration'],
           pal=colvec,breaks=breaks,
           main=main,reset=FALSE)
      } else {
        plot((sobj$polysf[datlist$poly.idx,] %>%
                mutate(concentration=datlist$cmat[it,]))['concentration'],
             pal=colvec,breaks=breaks,
             main=main,reset=FALSE)
      }
      mtext(ifelse(logscale,"Log scale","Linear scale"), side=3, line=0, adj=0, cex=0.8)
      if(fill.coast) plot(st_geometry(sobj$coastsf), col=coast.col, add=TRUE)
      if(make.movie) {
        box()
        if(date.axis) {
          par(usr=c(as.numeric(week.as.date(datlist$tvec[c(it1,it2)])),0,1))
          axis.Date(1)
          points(week.as.date(week),0,pch=16,cex=2,xpd=TRUE)
        }
        mtext(format(week.as.date(week), "%d-%m-%Y"), side=3, line=0, adj=1, cex=0.8)
        dev.off()
      }
    }
    if(make.movie) {
      image_list <- lapply(list.files(outdir,full.names=TRUE), image_read)
      poly_animated <- image_animate(image_join(image_list), fps=movie.fps)
      image_write(poly_animated, fmovie)
      if(movie.clean) {
        unlink(list.files(outdir,full.names=TRUE))
        unlink(outdir, recursive=TRUE)
      }
    }

  } else {
    warning(paste0("Plot type ",type," not recognised."))
  }
  invisible()
}


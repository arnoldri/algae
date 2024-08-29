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


#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' Test function - multiply by 2
//'
//' @description Just testing
//'
//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

//' @title Create mean matrix mumat
//'
//' @description Create the mean matrix mumat
//'
//' @export
// [[Rcpp::export]]
NumericMatrix rcpp_make_mumat(List parlist, List cpar) {

  int ne = parlist["ne"];
  int nt = cpar["nt"];
  int nm = cpar["nm"];

  NumericMatrix amat(nt+1,nm);
  NumericMatrix bmat(nt+1,nm);
  NumericMatrix mmat(nt+1,nm);
  NumericMatrix mumat(nt,nm);

  /*
  if(ne>0) {
    for(int k=0; k<ne; k++) {
       for(int i=0; i<wvec[k]; i++) {
         amat[tvec[k],]
       }
    }
  }
   */

  /*
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
  */

  return mumat;
}

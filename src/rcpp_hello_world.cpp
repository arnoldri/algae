
#include <Rcpp.h>
using namespace Rcpp;

//' Hello world function, this time in C
//'
//' @description Just getting going
//'
//' @export
// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z;
}

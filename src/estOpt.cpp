# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export]] 
List estOpt(NumericVector y,
              NumericVector r,
              NumericVector trt,
              NumericVector prtx,
              NumericVector rule,
              double tau,
              double lambda){
  
  double eff = 0;
  double risk = 0;
  double denom = 0;
  int ylength = y.length();
  int rlength = r.length();
  
  for(int i = 0; i < ylength; i++){
    eff += y[i] * (trt[i] * rule[i] + (1-trt[i]) * (1-rule[i])) / prtx[i];
    denom += (trt[i] * rule[i] + (1-trt[i]) * (1-rule[i])) / prtx[i];
  }
  for(int i = 0; i < rlength; i++){
    risk += r[i] * (trt[i] * rule[i] + (1-trt[i]) * (1-rule[i])) / prtx[i];
  }
  
  double output = eff / denom - lambda * (risk / denom - tau);
  return List::create(
    _["optValue"] = output,
    _["efficacy"] = eff / denom, 
    _["risk"] = risk / denom
  );
}

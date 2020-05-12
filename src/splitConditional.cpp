# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export]] 
Rcpp::List splitConditional(NumericVector zcut,
                            NumericMatrix zcutCat,
                            List datMatrix,
                            List parameters) {
  
  // Extract information from data matrix
  double useOtherNodes = as<double>(parameters["useOtherNodes"]);
  double isCtg = as<double>(parameters["isCtg"]);
  NumericVector y = as<NumericVector>(datMatrix["y"]);
  NumericVector x = as<NumericVector>(datMatrix["x"]);
  NumericVector ae = as<NumericVector>(datMatrix["ae"]);
  NumericVector prtx = as<NumericVector>(datMatrix["prtx"]);
  NumericVector trt = as<NumericVector>(datMatrix["trt"]);
  NumericVector trtNew = as<NumericVector>(datMatrix["trtNew"]);
  NumericVector inNode = as<NumericVector>(datMatrix["inNode"]);
  
  // Extract information from parameter list; 
  double nodeSize = as<double>(parameters["nodeSize"]);
  double trtSize = as<double>(parameters["trtSize"]);
  double maxScore = as<double>(parameters["maxScore"]);
  double maxRisk = as<double>(parameters["maxRisk"]);
  double lambda = as<double>(parameters["lambda"]);
  
  // Define several size parameters
  double totalSampSize = y.length();
  double zlength = zcutCat.ncol()*isCtg + zcut.length()*(1-isCtg);
  double ylength = y.length();
  double parentSize = 0;
  for(int i = 0; i < totalSampSize; i++){
    parentSize += inNode[i];
  }
  double divisor = totalSampSize*useOtherNodes + parentSize*(1-useOtherNodes);
  NumericVector outputL(zlength);
  NumericVector outputR(zlength);
  NumericVector outputLeft(zlength);
  NumericVector outputRight(zlength);
  NumericVector outputL2(zlength);
  NumericVector outputR2(zlength);
  
  // =========================================================
  // Candidate splits sending trt = 1 to the left node
  // =========================================================
  
  NumericVector condA(zlength);
  NumericVector condB1(zlength);
  NumericVector condB2(zlength);
  NumericVector condRiskL(zlength);
  NumericVector condRiskR(zlength);
  NumericVector condRiskL2(zlength);
  NumericVector condRiskR2(zlength);
  
  // Determine the pass status for each candidate split 
  // condA: number of obs per daughter node 
  // condB1: enough treated in left daughter
  // condB2: enough treated in right daughter
  // condRisk: Computes risk value function
  
  for(int i = 0; i < zlength; i++){
    for(int j = 0; j < ylength; j++){
      if(inNode[j]){
        if(isCtg){
          condA[i] += (zcutCat(j,i));
          condB1[i] += (zcutCat(j,i))*trt[j];
          condB2[i] += (1-zcutCat(j,i))*trt[j];
          condRiskL[i] += ae[j]*(trt[j]*(zcutCat(j,i)) / (divisor*prtx[j]) + 
            (1 - trt[j])*(1-zcutCat(j,i)) / (divisor*(1 - prtx[j])));
          condRiskR[i] += ae[j]*(trt[j]*(1-zcutCat(j,i)) / (divisor*prtx[j]) + 
            (1 - trt[j])*(zcutCat(j,i)) / (divisor*(1 - prtx[j])));
        } else{
          condA[i] += (x[j] <= zcut[i]);
          condB1[i] += (x[j] <= zcut[i])*trt[j];
          condB2[i] += (x[j] > zcut[i])*trt[j];
          condRiskL[i] += ae[j]*(((trt[j]*(x[j] <= zcut[i])) / prtx[j]) + (((1-trt[j])*(x[j] > zcut[i])) / (1-prtx[j])));
          condRiskL2[i] += (((trt[j]*(x[j] <= zcut[i])) / prtx[j]) + (((1-trt[j])*(x[j] > zcut[i])) / (1-prtx[j])));
          condRiskR[i] += ae[j]*(((trt[j]*(x[j] > zcut[i])) / prtx[j]) + (((1-trt[j])*(x[j] <= zcut[i])) / (1-prtx[j])));
	  condRiskR2[i] += (((trt[j]*(x[j] > zcut[i])) / prtx[j]) + (((1-trt[j])*(x[j] <= zcut[i])) / (1-prtx[j])));
        }
      } else{
        if(useOtherNodes){
          condRiskL[i] += ae[j]*(((trt[j]*trtNew[j]) / prtx[j]) + (((1-trt[j])*(1-trtNew[j])) / (1-prtx[j])));
          condRiskL2[i] += (((trt[j]*trtNew[j]) / prtx[j]) + (((1-trt[j])*(1-trtNew[j])) / (1-prtx[j])));
          condRiskR[i] += ae[j]*(((trt[j]*trtNew[j]) / prtx[j]) + (((1-trt[j])*(1-trtNew[j])) / (1-prtx[j])));
          condRiskR2[i] += (((trt[j]*trtNew[j]) / prtx[j]) + (((1-trt[j])*(1-trtNew[j])) / (1-prtx[j])));
        }
      }
    }
  }
  
  for(int i = 0; i < zlength; i++){
    if((condA[i] >= nodeSize) &&
       (parentSize - condA[i] >= nodeSize) &&
       (condB1[i] >= trtSize) &&
       (condB2[i] >= trtSize) &&
       (condA[i] - condB1[i] >= trtSize) &&
       (parentSize - condA[i] - condB2[i] >= trtSize)){
      // if(condRiskL[i] <= maxRisk*condRiskL2[i]){
      for(int j = 0; j < ylength; j++){
        if(inNode[j]){
          if(isCtg){
            outputL[i] += y[j]*(trt[j]*(zcutCat(j,i)) / (prtx[j]) + (1 - trt[j])*(1-zcutCat(j,i)) / (1 - prtx[j]));
	    outputL2[i] += (trt[j]*(zcutCat(j,i)) / (divisor*prtx[j]) + (1 - trt[j])*(1-zcutCat(j,i)) / (divisor*(1 - prtx[j])));
          } else{
            outputL[i] += (y[j])*((trt[j]*(x[j] <= zcut[i]) / prtx[j]) + ((1-trt[j])*(x[j] > zcut[i]) / (1-prtx[j])));
	    outputL2[i] += ((trt[j]*(x[j] <= zcut[i]) / prtx[j]) + ((1-trt[j])*(x[j] > zcut[i]) / (1-prtx[j])));
          }
        } else{
          if(useOtherNodes){
            outputL[i] += (y[j])*(((trt[j]*trtNew[j]) / prtx[j]) + (((1-trt[j])*(1-trtNew[j])) / (1-prtx[j])));
	    outputL2[i] += (((trt[j]*trtNew[j]) / prtx[j]) + (((1-trt[j])*(1-trtNew[j])) / (1-prtx[j])));
          }
        }
      }
      // } else{
      //   outputL[i] = -1E20;
      // }
      
      // if(condRiskR[i] <= maxRisk*condRiskR2[i]){
      for(int j = 0; j < ylength; j++){
        if(inNode[j]){
          if(isCtg){
            outputR[i] += y[j]*(trt[j]*(1-zcutCat(j,i)) / (prtx[j]) + (1 - trt[j])*(zcutCat(j,i)) / (1 - prtx[j]));
	    outputR2[i] += (trt[j]*(1-zcutCat(j,i)) / (prtx[j]) + (1 - trt[j])*(zcutCat(j,i)) / (1 - prtx[j]));
          } else{
            outputR[i] += (y[j])*((trt[j]*(x[j] > zcut[i]) / prtx[j]) + (1-trt[j])*(x[j] <= zcut[i]) / (1-prtx[j]));
	    outputR2[i] += ((trt[j]*(x[j] > zcut[i]) / prtx[j]) + (1-trt[j])*(x[j] <= zcut[i]) / (1-prtx[j]));
          }
        } else{
          if(useOtherNodes){
            outputR[i] += (y[j])*((trt[j]*trtNew[j] / prtx[j]) + (1-trt[j])*(1-trtNew[j]) / (1-prtx[j]));
	    outputR2[i] += ((trt[j]*trtNew[j] / prtx[j]) + (1-trt[j])*(1-trtNew[j]) / (1-prtx[j]));
          }
        }
      }
      // } else{
      //   outputR[i] = -1E20;
      // }
    }
  }
  
  for(int i = 0; i < zlength; i++){
    outputLeft[i] = (outputL[i] / outputL2[i]) - lambda * ((condRiskL[i] / condRiskL2[i]) - maxRisk);
    outputRight[i] = (outputR[i] / outputR2[i]) - lambda * ((condRiskR[i] / condRiskR2[i]) - maxRisk);
  }
  
  // Get output direction for splits
  Rcpp::StringVector outputDirection(zlength);
  NumericVector output(zlength);
  for(int i = 0; i < zlength; i++){
    if((outputLeft[i] > outputRight[i]) && (outputLeft[i] >= maxScore)){
      outputDirection[i] = "l";
      output[i] = outputLeft[i];
    } else if((outputLeft[i] < outputRight[i]) && (outputRight[i] >= maxScore)){
      outputDirection[i] = "r";
      output[i] = outputRight[i];
    } else{
      outputDirection[i] = "neither";
      output[i] = -1E10;
    }
  }
  
  // Compile output
  return List::create(
    _["direction"] = outputDirection,
    _["output"] = output, 
    _["riskL"] = condRiskL / condRiskL2, 
    _["riskR"] = condRiskR / condRiskR2, 
    _["valueL"] = outputL / outputL2, 
    _["valueR"] = outputR / outputR2, 
    _["divisor"] = divisor,
    _["zcut"] = zcut
  );
}

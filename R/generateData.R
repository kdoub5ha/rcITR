#' @title Example data for rcITR package demonstration
#' @description Generates data for demonstration of rcDT and rcRF model construction.
#' @param n numeric. Number of example samples to be generated. Defaults to 1000.
#' @return data.frame with 10 covariates (x1-x10), efficacy (y), risk (r), 
#'         treatment indicator (trt), and propensity score (prtx).
#' @import randomForest
#' @export
#' @examples
#' dat <- generateData()

generateData <- function(n = 1000){
  trt <- sample(c(-1,1), n, replace = TRUE)
  X <- matrix(round(runif(10*n), 3), n, 10, dimnames = list(c(1:n), paste0("x", 1:10)))
  subgroup <- (X[,1] <= 0.7 & X[,2] <= 0.7)
  Y <- 1 - X[,3] + X[,4] + 3*subgroup*(trt == 1) + 3*(1-subgroup)*(trt == -1) + rnorm(n)
  R <- 2 + X[,4] + 2*X[,3] + (2*X[,2] + 0.1)*trt + rnorm(n, sd = 0.25)
  out <- data.frame(X, 
                    trt = ifelse(trt == -1, 0, 1), 
                    y = round(Y, 4), 
                    r = round(R, 4), 
                    prtx = 0.5)
  return(out)
}
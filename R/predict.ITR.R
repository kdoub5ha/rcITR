#' @title Treatment prediction for rcDT and rcRF models
#'
#' @description Returns treatment predictions for a rcDT and rcRF models given a new data set. If the 
#' input is rcRF (forest), then the proportion of trees voting for treatment (`trt=1`) is returned. 
#' If the input is rcDT (single tree), then the function returns the vote (0 / 1) for the model. 
#' 
#' @param fit tree or forest object from `rcDT` or `rcRF`.
#' @param new.data data.frame of new observations
#' @param split.var numeric vector indicating columns of covariates
#' @param ctgs numeric vector of columns of categorical covariates. 
#' @return A list of prediction summaries
#' @return \item{SummaryTreat}{proportion of trees voting for treatment (trt=1). 
#' If input is rcDT (single tree) then SummaryTreat is a single number. 
#' If input is rcRF (forest) then SummaryTreat is a vector equal to the length of the number of trees.}
#' @return \item{trt.pred}{vector of treatment assignments {0, 1} based on the tree vote (single tree) or majority of tree votes (forest). This vector has length equal to the number of rows in `new.data`.}
#' @return \item{n.trees}{number of tree in `fit`}
#' @return \item{tree.votes}{matrix of votes for each tree for each subject in `new.data`. Rows correspond to trees in `fit` and columns correspond to subjects in `new.dat`.}
#' @return \item{data}{input data frame `new.data`}
#' @return \item{NA.trees}{number of trees returning no votes. In a forest, this is the number of null trees.}
#' @export
#' @examples
#' # Generate simulated data
#' set.seed(123)
#' dat <- generateData()
#' 
#' # Generates rcDT using simualated data with splitting variables located in columns 1-10.
#' rcDT.fit <- rcDT(data = dat, 
#'                  split.var = 1:10, 
#'                  risk.threshold = 2.75, 
#'                  lambda = 1)
#' # Predict treatment assignments for 1000 observations in `dat` using the rcDT model
#' preds.rcDT <- predict.ITR(fit = rcDT.fit, new.data = dat, split.var = 1:10)
#' 
#' # Generates rcRF using simualated data with splitting variables located in columns 1-10.
#' set.seed(2)
#' rcRF.fit <- rcRF(data = dat, 
#'                  split.var = 1:10, 
#'                  ntree = 200,
#'                  risk.threshold = 2.75, 
#'                  lambda = 1)
#' # Predict treatment assignments for 1000 observations in `dat` using the rcRF model
#' preds.rcRF <- predict.ITR(fit = rcRF.fit, new.data = dat, split.var = 1:10)
#' 

predict.ITR <- function(fit, 
                        new.data,  
                        split.var, 
                        ctgs = NULL){
  if(is.null(dim(fit))){
    trees <- fit$TREES
    n.trees <- length(trees)
  } else{
    trees <- fit
    n.trees <- 1
  }
  dat <- new.data
  n <- nrow(dat)
  out <- NULL
 
  result <- sapply(1:n.trees, function(i){
    if(is.null(dim(fit))){
      tre <- trees[[i]]
    } else{
      tre <- trees
    }
    
    if(nrow(tre) > 0){
      if(!is.na(tre[1,6])){
        idx <- !is.na(tre$cut.2)
        cutPoint <- as.numeric(tre$cut.2[idx])
        splitVar <- as.integer(tre$var[idx])
        treNodes <- as.character(tre$node[idx])
        direction <- as.character(tre$cut.1[idx])
        Data <- as.matrix(dat[,split.var,drop=F])
        send <- SendDown(cutPoint, splitVar, Data, treNodes, direction)
        trt.pred <- send$trt.pred
      } else{
        trt.pred <- rep(NA, n)
      }
    } else{
      trt.pred <- rep(NA, n)
    }
    
    return(trt.pred) 
  })
  
  if(!(is.null(dim(result)) & length(result) == 1)){
    out$SummaryTreat <- apply(result, 1, FUN = mean, na.rm=T)
  } else{
    out$SummaryTreat <- result
  }
  if(is.null(dim(fit))){
    out$trt.pred <- ifelse(out$SummaryTreat < 0.5, 0, 1)
  } else{
    out$trt.pred <- out$SummaryTreat
  }

  out$n.trees <- n.trees
  out$tree.votes <- result
  out$data <- new.data
  out$NA.trees <- ifelse((!(is.null(dim(result)) & length(result) == 1)), sum(is.na(result[1,])), NA)
  return(out)
}
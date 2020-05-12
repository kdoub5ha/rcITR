#' @title Builds an rcRF model for risk controlled optimization
#' 
#' @description This function constructs a random forest of rcDT trees method. 
#' A forest object can be an argument into `predict.ITR()` along with data in order to 
#' obtain treatment predictions. An output from this function can also be given to the `Variable.Importance.ITR()` 
#' function to estimate predictor importance. 
#' 
#' @param data data.frame. Data used to construct rcRF model.  
#' Must contain efficacy variable (y), 
#' risk variable (r), 
#' binary treatment indicator coded as 0 / 1 (trt), 
#' propensity score (prtx),
#' candidate splitting covariates.
#' @param split.var numeric vector. Columns of spliting variables.
#' @param efficacy char. Efficacy outcome column. Defaults to 'y'.
#' @param risk char. Risk outcome column. Defaults to 'r'.
#' @param col.trt char. Treatment column name
#' @param col.ptrx char. Propensity score column name.
#' @param risk.control logical. Should risk be controlled? Defaults to TRUE.
#' @param risk.threshold numeric. Desired level of risk control. 
#' @param lambda numeric. Penalty parameter for risk scores. Defaults to 0, i.e. no constraint.
#' 
#' Optional arguments
#' @param test data.frame of testing observations. Should be formatted the same as 'data'.
#' @param N0 numeric specifying minimum number of observations required to call a node terminal. Defaults to 20.
#' @param n0 numeric specifying minimum number of treatment/control observations needed in a split to declare a node terminal. Defaults to 5. 
#' @param max.depth numeric specifying maximum depth of the tree. Defaults to 15 levels. 
#' @param mtry numeric specifying the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @param ntree numeric. Number of trees generated. Defaults to 500.
#' @param stabilize logical indicating if efficacy should be modeled using residuals. Defaults to TRUE. 
#' @param stabilize.type character specifying method used for estimating residuals. Current options are 'linear' for linear model (default) and 'rf' for random forest. 
#' @param use.other.nodes logical. Should global estimator of objective function be used. Defaults to TRUE. 
#' @param ctg numeric vector corresponding to the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param avoid.nul.tree logical. Should null trees be discarded?
#' @param verbose logical. Give updates about forest progression?
#' @param AIPWE logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. Not available yet. 
#' @param extremeRandomized logical. Experimental for randomly selecting cutpoints in a random forest model. Defaults to FALSE and users should change this at their own peril. 
#' @param order.importances logical. Should importances be ordered (if requested)?
#' @param importance.measures logical. Indicated if variable importance measures should be estimated and returned. Defaults to FALSE.
#' @return A list of characteristics of the forest.
#' @return \item{ID.Boots.Samples}{list of bootstrap sample IDs}
#' @return \item{TREES}{list of trees}
#' @return \item{Model.Specification}{information about the input parameters of the forest}
#' @return \item{...}{Summaries for in and out of bag samples}
#' @import randomForest
#' @export
#' @examples
#' set.seed(123)
#' dat <- generateData(n = 500)
#' # Generates rcRF model using simualated data with splitting variables located in columns 1-10.
#' fit <- rcRF(dat = dat, split.var = 1:10, ntree = 200,
#'             risk.control = TRUE, risk.threshold = 2.75, 
#'             lambda = 1)


rcRF <- function(data, 
                 split.var, 
                 efficacy = "y",
                 risk = "r",
                 col.trt = "trt", 
                 col.prtx = "prtx", 
                 risk.control = TRUE, 
                 risk.threshold = NA, 
                 lambda = 0,
                 stabilize = TRUE, 
                 stabilize.type = c('linear', 'rf'), 
                 test = NULL, 
                 ctg = NULL,
                 N0 = 20, 
                 n0 = 5,  
                 max.depth = 10,
                 ntree = 500, 
                 mtry = max(floor(length(split.var)/3), 1),
                 avoid.nul.tree = FALSE, 
                 AIPWE = FALSE, 
                 verbose = FALSE, 
                 use.other.nodes = TRUE, 
                 extremeRandomized = FALSE, 
                 importance.measures = FALSE, 
                 order.importances = TRUE)
{
  
  require(randomForest)
  
  # input checks
  if(!is.data.frame(data)) stop("data argument must be dataframe")
  if(!is.numeric(split.var)) stop("split.var must be numeric vector")
  
  if(!is.character(efficacy)) stop("efficacy argument must be character")
  if(!efficacy %in% colnames(data)) stop("efficacy argument is not in data")
  if(!is.character(risk)) stop("risk argument must be character")
  if(!risk %in% colnames(data)) stop("risk argument is not in data")
  
  if(risk.control & is.na(risk.threshold)) warning("risk.contrl is TRUE, but risk.threshold not specified")
  if(!any(stabilize.type %in% c("linear", "rf"))) stop("linear and rf values supported for stabilize.type")
  if(mtry > length(split.var)){
    warning("mtry is larger than split.var length -- setting mtry to length(split.var)")
    mtry <- length(split.var)
  }
  
  stabilize.type <- match.arg(stabilize.type)
  
  if(sum(data$trt %in% c(0,1)) != nrow(data)){
    data$trt <- ifelse(data$trt == -1, 0, 1)
    message("Assuming trt indicator is of form -1/+1 and changed values to 0/1")
  }
  
  out <- as.list(NULL)
  out$ID.Boots.Samples  <- as.list(1:ntree)
  out$TREES <- as.list(1:ntree)
  out$preds.oob <- matrix(NA, nrow = nrow(data), ncol = ntree)
  out$preds.inbag <- matrix(NA, nrow = nrow(data), ncol = ntree)
  out$preds.cumulative.oob <- matrix(NA, nrow = nrow(data), ncol = ntree)
  out$preds.cumulative.inbag <- matrix(NA, nrow = nrow(data), ncol = ntree)
  out$value.oob <- rep(NA, ntree)
  out$value.inbag <- rep(NA, ntree)
  out$risk.oob <- rep(NA, ntree)
  out$risk.inbag <- rep(NA, ntree)
  out$risk.bound.values <- NULL
  
  # Initialize output if test set is included
  if(!is.null(test)) {
    out$test.preds <- matrix(NA, nrow = nrow(test), ncol = ntree)
    out$test.preds.cumulative <- matrix(NA, nrow = nrow(test), ncol = ntree)
    out$test.value <- rep(NA, ntree)
    out$test.risk <- rep(NA, ntree)
  }
  
  # set inputs
  stabilize.type <- match.arg(stabilize.type)
  if(risk.control){
    if(is.na(risk.threshold)) stop("Risk allowance level must be specified numeric")
  }
  
  # set parameters for splitting criteria
  b <- 1
  while(b <= ntree){
    # TAKE BOOTSTRAP SAMPLES
    id.b <- sample(1:nrow(data), size=nrow(data), replace = TRUE)
    dat.b <- data[id.b,]
    dat.test <- data[-unique(id.b),]
    
    # Generate tree based on b-th bootstrap sample
    tre.b <- rcDT(data = dat.b, 
                  split.var = split.var, 
                  test = test, 
                  min.ndsz = N0, 
                  n0 = n0, 
                  efficacy = efficacy, 
                  risk = risk, 
                  col.trt = col.trt,
                  col.prtx = col.prtx,
                  lambda = lambda,
                  risk.control = risk.control, 
                  risk.threshold = risk.threshold, 
                  stabilize = stabilize,
                  stabilize.type = stabilize.type,
                  ctg = ctg, 
                  max.depth = max.depth, 
                  AIPWE = AIPWE, 
                  mtry = mtry,
                  use.other.nodes = use.other.nodes, 
                  extremeRandomized = extremeRandomized)
    
    if(avoid.nul.tree) {
      if(nrow(tre.b$tree) > 1) {
        out$ID.Boots.Samples[[b]] <- id.b
        out$TREES[[b]] <- tre.b$tree
        b <- b + 1
      }
    } else {
      out$ID.Boots.Samples[[b]] <- id.b
      out$TREES[[b]] <- tre.b$tree
      if(nrow(tre.b$tree) > 1){
        preds <- predict.ITR(tre.b$tree, data, split.var)$trt.pred
        inbag.idx <- seq_along(1:nrow(data)) %in% unique(id.b)
        out$preds.oob[,b] <- ifelse(inbag.idx, NA, preds)
        out$preds.inbag[,b] <- ifelse(inbag.idx, preds, NA)
      }
      out$preds.cumulative.oob[,b] <- ifelse(rowMeans(out$preds.oob[,1:b,drop=F], na.rm = T) > 0.5, 1, 0)
      out$preds.cumulative.inbag[,b] <- ifelse(rowMeans(out$preds.inbag[,1:b,drop=F], na.rm = T) > 0.5, 1, 0)
      
      if(b > 5){
        out$value.inbag[b] <- sum(data[,efficacy] * (out$preds.cumulative.inbag[,b] == data[,col.trt]) / data[,col.prtx], na.rm = T) / 
          sum((out$preds.cumulative.inbag[,b] == data[,col.trt]) / data[,col.prtx], na.rm = T)
        out$value.oob[b] <- sum(data[,efficacy] * (out$preds.cumulative.oob[,b] == data[,col.trt]) / data[,col.prtx], na.rm = T) / 
          sum((out$preds.cumulative.oob[,b] == data[,col.trt]) / data[,col.prtx], na.rm = T)
        out$risk.inbag[b] <- sum(data[,risk] * (out$preds.cumulative.inbag[,b] == data[,col.trt]) / data[,col.prtx], na.rm = T) / 
          sum((out$preds.cumulative.inbag[,b] == data[,col.trt]) / data[,col.prtx], na.rm = T)
        out$risk.oob[b] <- sum(data[,risk] * (out$preds.cumulative.oob[,b] == data[,col.trt]) / data[,col.prtx], na.rm = T) / 
          sum((out$preds.cumulative.oob[,b] == data[,col.trt]) / data[,col.prtx], na.rm = T)
      }
      
      if(!is.null(test)){
        if(nrow(tre.b$tree) > 1){
          out$test.preds[,b] <- predict.ITR(tre.b$tree, as.data.frame(test), split.var)$trt.pred
        }
        if(b > 5){
          out$test.preds.cumulative[,b] <- ifelse(rowMeans(out$test.preds[,1:b,drop=F], na.rm = T) > 0.5, 1, 0)
          out$test.value[b] <- sum(test[,efficacy] * (out$test.preds.cumulative[,b] == test[,col.trt]) / test[,col.prtx], na.rm = T) / 
            sum((out$test.preds.cumulative[,b] == test[,col.trt]) / test[,col.prtx], na.rm = T)
          out$test.risk[b] <- sum(test[risk] * (out$test.preds.cumulative[,b] == test[,col.trt]) / test[,col.prtx], na.rm = T) / 
            sum((out$test.preds.cumulative[,b] == test[,col.trt]) / test[,col.prtx], na.rm = T)
        }
      }
      
      b <- b + 1
    }
    if(verbose){
      if(b %in% seq(0, ntree, 100)) print(paste0("On Tree ", b))
    }
  }
  
  Model.Specification <- as.list(NULL)
  Model.Specification$data <- data
  Model.Specification$split.var <- split.var
  Model.Specification$ctg <- ctg
  Model.Specification$efficacy <- efficacy
  Model.Specification$risk <- risk
  Model.Specification$col.trt <- col.trt
  Model.Specification$col.prtx <- col.prtx
  Model.Specification$lambda <- lambda
  Model.Specification$risk.threshold <- risk.threshold
  out$Model.Specification <- Model.Specification
  
  if(importance.measures){
    importances <- Variable.Importance.ITR(out, sort = order.importances)
    out$importances <- importances
  }
  return(out)
}
#' @title rcRF model selection
#' @description Performs model selection for rcRF model to 
#' select the best penalty parameter ($\lambda$).
#' @param data data.frame. Data used to construct rcDT model.  
#' Must contain efficacy variable (y), 
#' risk variable (r), 
#' binary treatment indicator coded as 0 / 1 (trt), 
#' propensity score (prtx),
#' candidate splitting covariates.
#' @param split.var numeric vector. Columns of spliting variables.
#' @param efficacy char. Efficacy outcome column. Defaults to 'y'.
#' @param risk char. Risk outcome column. Defaults to 'r'.
#' @param col.trt char. Treatment indicator column name. Should be of form 0/1 or -1/+1.
#' @param col.prtx char. Propensity score column name. 
#' @param risk.control logical. Should risk be controlled? Defaults to TRUE.
#' @param risk.threshold numeric. Desired level of risk control. 
#' @param ntree numeric. Number of trees to construct.
#' @param test data.frame of testing observations. Should be formatted the same as 'data'.
#' @param N0 numeric specifying minimum number of observations required to call a node terminal. Defaults to 20.
#' @param n0 numeric specifying minimum number of treatment/control observations needed in a split to declare a node terminal. Defaults to 5. 
#' @param max.depth numeric specifying maximum depth of the tree. Defaults to 15 levels. 
#' @param mtry numeric specifying the number of randomly selected splitting variables to be included. 
#' Defaults to the greater of 1 and length(split.var)/3.
#' @param stabilize logical indicating if efficacy should be modeled using residuals. Defaults to TRUE. 
#' @param stabilize.type character specifying method used for estimating residuals. Current options are 'linear' for linear model (default) and 'rf' for random forest. 
#' @param use.other.nodes logical. Should global estimator of objective function be used. Defaults to TRUE. 
#' @param ctg numeric vector corresponding to the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param AIPWE logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. Not available yet. 
#' @param extremeRandomized logical. Experimental for randomly selecting cutpoints in a random forest model. Defaults to FALSE and users should change this at their own peril. 
#' @param lambda.upper numeric. Upper bound for risk penalty. An attempt at reasonable selection will be performed automatically. 
#' @param importance logical. Indicated if variable importance measures should be estimated and returned. Defaults to FALSE.
#' @param order.importances logical. Should importances be ordered (if requested)?
#' @param max.iter numeric. Indicates the maximum number of forest iterations to perform. Defaults to 10.
#' @param risk.tolerance numeric. Two component vector giving the bound on risk that is acceptable (acceptable risk range is calcuated as risk.threshold * risk.tolerance). Defaults to c(0.995, 1.005), i.e. 0.5\% tolerance. 
#' @param avoid.nul.tree logical. Should null trees be discarded?
#' @param verbose logical. Give updates about forest progression?
#' @return A summary of the cross validation including optimal penalty parameter and the optimal model. 
#' @return \item{best.fit}{optimal rcRF model}
#' @return \item{lambda}{optimal lambda value selected}
#' @return \item{oob.risk}{out-of-bag risk from best model}
#' @return \item{converged}{max number of iterations reached?}
#' @return \item{importances}{importance measures, if requested}
#' @return \item{risks}{vector of risk scores obtained over tuning procedure}
#' @return \item{lambdas}{vector of lambda values tried over tuning procedure}
#' @return \item{time.elapsed}{elapsed time for model tuning}
#' @import randomForest
#' @export
#' @examples
#' 
#' # Grow large tree
#' set.seed(123)
#' dat <- generateData()
#' fit <- rcRF.select(data = dat, 
#'                    split.var = 1:10,
#'                    risk.threshold = 2.75)
#' 


rcRF.select <- function(data, 
                        split.var,
                        test = NULL,
                        N0 = 20, 
                        n0 = 5, 
                        efficacy = 'y',
                        risk = 'r',
                        col.trt = "trt",
                        col.prtx = "prtx",
                        ntree = 500,
                        lambda.upper = NA,
                        risk.control = TRUE, 
                        risk.threshold = NA, 
                        AIPWE = FALSE, 
                        ctg = NA, 
                        mtry = max(floor(length(split.var)/3), 1),
                        avoid.nul.tree = FALSE, 
                        max.depth = 15,
                        stabilize.type = c('linear', 'rf'), 
                        stabilize = TRUE, 
                        verbose = FALSE, 
                        use.other.nodes = TRUE, 
                        extremeRandomized = FALSE, 
                        importance = FALSE, 
                        order.importances = TRUE,
                        max.iter = 10,
                        risk.tolerance = c(0.995, 1.005)){

  start.time <- Sys.time()
  # input checks
  if(!is.data.frame(data)) stop("data argument must be dataframe")
  if(!is.numeric(split.var)) stop("split.var must be numeric vector")
  
  if(!is.character(efficacy)) stop("efficacy argument must be character")
  if(!efficacy %in% colnames(data)) stop("efficacy argument is not in data")
  if(!is.character(risk)) stop("risk argument must be character")
  if(!risk %in% colnames(data)) stop("risk argument is not in data")
  if(!is.numeric(max.iter) | length(max.iter) != 1) stop("max.iter must be numeric of length 1")
  if(length(risk.tolerance) != 2 | !is.numeric(risk.tolerance)) 
    stop("risk.tolerance must be two component numeric vector")
  
  if(risk.control & is.na(risk.threshold)) warning("risk.contrl is TRUE, but risk.threshold not specified")
  if(!any(stabilize.type %in% c("linear", "rf"))) stop("linear and rf values supported for stabilize.type")
  if(mtry > length(split.var)){
    warning("mtry is larger than split.var length -- setting mtry to length(split.var)")
    mtry <- length(split.var)
  }
  
  stabilize.type <- match.arg(stabilize.type)
  
  if(sum(data[,col.trt] %in% c(0,1)) != nrow(data)){
    data[,col.trt] <- ifelse(data[,col.trt] == -1, 0, 1)
    message("Assuming trt indicator is of form -1/+1 and changed values to 0/1")
  }
  
  # retrieve lambda information (if provided)
  if(is.na(lambda.upper)){
    # define range of risk penalty values
    EY.0 <- mean(data[data[,col.trt] == 0,efficacy] - mean(data[,efficacy])) 
    EY.1 <- mean(data[data[,col.trt] == 1,efficacy] - mean(data[,efficacy])) 
    ER.0 <- mean(data[data[,col.trt] == 0,risk]) 
    ER.1 <- mean(data[data[,col.trt] == 1,risk]) 
    
    # define upper bound for lambda
    lambda.upper <- 1.25*(EY.1 - EY.0 + ER.0) / (ER.1 - risk.threshold)
  } else{
    if(!is.numeric(lambda.upper)) stop("lambda upper bound must be numeric")
  }
  
  fits <- vector("list")
  
  # scan over lambda values to select final lambda value
  lambdas <- 0.95*lambda.upper
  iter <- 1
  risks <- max(
    c(sum(data[,risk] * (data[,col.trt] == 0) / data[,col.prtx]) / sum((data[,col.trt] == 0) / data[,col.prtx]),
      sum(data[,risk] * (data[,col.trt] == 1) / data[,col.prtx]) / sum((data[,col.trt] == 1) / data[,col.prtx]))
  )

  while(((risks[length(risks)] > risk.tolerance[2]*risk.threshold) | 
         (risks[length(risks)] < risk.tolerance[1]*risk.threshold))& 
        iter < max.iter + 1){
  
    fit.lam <- rcRF(data = data, 
                    split.var = split.var, 
                    efficacy = efficacy, 
                    risk = risk, 
                    col.trt = col.trt, 
                    col.prtx = col.prtx, 
                    risk.control = risk.control, 
                    risk.threshold = risk.threshold, 
                    lambda = lambdas[length(lambdas)],
                    stabilize = stabilize,
                    stabilize.type = stabilize.type, 
                    test = test, 
                    ctg = ctg,
                    N0 = N0, 
                    n0 = n0, 
                    max.depth = max.depth,
                    ntree = ntree, 
                    mtry = mtry, 
                    avoid.nul.tree = avoid.nul.tree, 
                    AIPWE = AIPWE, 
                    verbose = FALSE,
                    use.other.nodes = use.other.nodes, 
                    extremeRandomized = extremeRandomized, 
                    importance = FALSE, 
                    order.importances = order.importances)
    fits <- c(fits, list(fit.lam))
    risks <- c(risks, fit.lam$risk.oob[length(fit.lam$risk.oob)])
    if(verbose){
      cat("Completed Iteration", iter, "; OOB risk =", risks[length(risks)], "\n")
    }
    
    iter <- iter + 1

    if(length(lambdas) == 1){
      if(risks[length(risks)] > risk.threshold){
        lambdas <- c(lambdas, 1.5*lambdas[length(lambdas)])
      } else{
        lambdas <- c(lambdas, 0.5*lambdas[length(lambdas)])
      }
    } else{
      if(risks[length(risks)] > risk.threshold){
        # current risk is too high, need to increase lambda
        if(length(which(risks[-1] < risks[length(risks)])) > 0){
          idx <- which(risks[-1] < risks[length(risks)])
          base.risk <- risks[-1][idx][which.min(abs(risks[-1][idx] - risks[length(risks)]))]
          idxx <- which(risks == base.risk)
          base.lambda <- lambdas[idxx-1]
          base.prop <- 1 - abs(risks[length(risks)] - risk.threshold) / abs(base.risk - risks[length(risks)])
          lambdas <- c(lambdas, base.lambda - base.prop*abs(base.lambda - lambdas[length(lambdas)]))
        } else{
          lambdas <- c(lambdas, 1.5 * lambdas[length(lambdas)]) 
        }
      } else{
        if(length(which(risks[-1] > risks[length(risks)])) > 0){
          idx <- which(risks[-1] > risks[length(risks)])
          base.risk <- risks[-1][idx][which.min(abs(risks[-1][idx] - risks[length(risks)]))]
          idxx <- which(risks == base.risk)
          base.lambda <- lambdas[idxx-1]
          base.prop <- 1 - abs(risks[length(risks)] - risk.threshold) / abs(base.risk - risks[length(risks)])
          lambdas <- c(lambdas, base.lambda + base.prop*abs(base.lambda - lambdas[length(lambdas)]))
        } else{
          lambdas <- c(lambdas, 0.5 * lambdas[length(lambdas)]) 
        }
      }
    }
  }

  efficacies <- do.call(c, lapply(fits, function(ff) ff$value.oob[length(ff$value.oob)]))
  risks <- risks[-1]
  best.idx1 <- which(risks < risk.threshold * risk.tolerance[2])
  best.idx2 <- which.max(efficacies[best.idx1])
  best.fit <- fits[best.idx1][[best.idx2]]
  best.lambda <- lambdas[best.idx1][best.idx2]
  best.risk <- risks[best.idx1][best.idx2]
  end.time <- Sys.time()

  if(importance){
    if(verbose) cat("Calculating Importances")
    importances <- Variable.Importance.ITR(rcRF.fit = best.fit, sort = order.importances)
  } else{
    importances <- NA
  }
  setNames(list(best.fit,
                best.lambda,
                best.risk,
                as.logical(iter <= max.iter),
                importances,
                risks[-1],
                lambdas[-length(lambdas)],
                end.time - start.time), 
           c("best.fit", 
             "lambda",
             "oob.risk",
             "converged",
             "importances",
             "risks",
             "lambdas",
             "elapsed.time"))
}

#' @title Optimal rcDT model selection
#' @description Performs k-fold cross validation for tuning of risk and tree size 
#' paramters to select the optimal rcDT model.
#' @param data data.frame. Data used to construct rcDT model.  
#' Must contain efficacy variable (y), 
#' risk variable (r), 
#' binary treatment indicator coded as 0 / 1 (trt), 
#' propensity score (prtx),
#' candidate splitting covariates (split.var).
#' @param split.var numeric vector. Columns of spliting variables.
#' @param efficacy char. Efficacy outcome column. Defaults to 'y'.
#' @param risk char. Risk outcome column. Defaults to 'r'.
#' @param col.trt char. Treatment indicator column name. Should be of form 0/1 or -1/+1.
#' @param col.prtx char. Propensity score column name. 
#' @param risk.control logical. Should risk be controlled? Defaults to TRUE.
#' @param risk.threshold numeric. Desired level of risk control. 
#' @param lambda.seq numeric vector. Identifies sequence of risk penalty parameters to be considered. 
#' Defaults to NA and will attempt to identify reasonable range. 
#' @param lambda.length numeric indicating number of risk penalty parameters to use in tuning. 
#' Larger values will cause model selection to be slower. Defaults to 50.
#' @param n.folds numeric. Number of folds to use in k-fold cross validation. Defaults to 10.
#' @param test data.frame of testing observations. Should be formatted the same as 'data'.
#' @param N0 numeric specifying minimum number of observations required to call a node terminal. Defaults to 20.
#' @param n0 numeric specifying minimum number of treatment/control observations needed in a split to declare a node terminal. Defaults to 5. 
#' @param max.depth numeric specifying maximum depth of the tree. Defaults to 15 levels. 
#' @param mtry numeric specifying the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @param stabilize logical indicating if efficacy should be modeled using residuals. Defaults to TRUE. 
#' @param stabilize.type character specifying method used for estimating residuals. Current options are 'linear' for linear model (default) and 'rf' for random forest. 
#' @param sort internal use.
#' @param use.other.nodes logical. Should global estimator of objective function be used. Defaults to TRUE. 
#' @param use.bootstrap logical. Should a bootstrap resampling be done? Defaults to FALSE.
#' @param ctg numeric vector corresponding to the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param AIPWE logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. Not available yet. 
#' @param extremeRandomized logical. Experimental for randomly selecting cutpoints in a random forest model. Defaults to FALSE and users should change this at their own peril. 
#' @param verbose logical. Should tuning progress bar be displayed. Defaults to TRUE.
#' @return A summary of the cross validation including optimal penalty parameter and the optimal model. 
#' @return \item{best.tree}{optimal rcDT model}
#' @return \item{alpha}{tree size penalty}
#' @return \item{lambda}{risk penalty}
#' @return \item{full.tree}{unpruned tree}
#' @return \item{pruned.tree}{output from pruning of `full.tree`}
#' @return \item{subtrees}{sequence of optimally pruned subtrees}
#' @return \item{best.tree.summaries}{summary across trees}
#' @return \item{in.train}{training samples from splits}
#' @return \item{in.test}{testing samples from splits}
#' @return \item{elapsed.time}{time elapsed during model tuning}
#' @import randomForest
#' @import pbapply
#' @export
#' @examples
#' 
#' # Grow large tree
#' set.seed(123)
#' dat <- generateData()
#' fit <- rcDT.select(data = dat, 
#'                    split.var = 1:10, 
#'                    nfolds = 5,
#'                    risk.threshold = 2.75)
#' 


rcDT.select <- function(data, 
                        split.var, 
                        N0 = 20, 
                        n0 = 5, 
                        efficacy = 'y',
                        risk = 'r',
                        col.trt = "trt",
                        col.prtx = "prtx",
                        lambda.seq = NA,
                        lambda.length = 50,
                        risk.control = TRUE, 
                        risk.threshold = NA, 
                        nfolds = 10, 
                        AIPWE = FALSE, 
                        sort = TRUE, 
                        ctg = NA, 
                        mtry = length(split.var),
                        max.depth = 15,
                        stabilize.type = c('linear', 'rf'), 
                        stabilize = TRUE, 
                        use.other.nodes = TRUE, 
                        use.bootstrap = FALSE,
                        extremeRandomized = FALSE,
                        verbose = TRUE){
  
  start.time <- Sys.time()
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

  # retrieve lambda information (if provided)
  if(is.na(lambda.seq)){
    # define range of risk penalty values
    EY.0 <- mean(data[data$trt == 0,efficacy] - mean(data[,efficacy])) 
    EY.1 <- mean(data[data$trt == 1,efficacy] - mean(data[,efficacy])) 
    ER.0 <- mean(data[data$trt == 0,risk]) 
    ER.1 <- mean(data[data$trt == 1,risk]) 
    
    # define upper bound for lambda
    lambda.upper <- 1.25*(EY.1 - EY.0 + ER.0) / (ER.1 - risk.threshold)
    lambdas <- round(seq(0, lambda.upper, length.out = lambda.length+1), 5)[-1]
  } else{
    if(length(lambda.seq) > 100){
      lambdas <- sort(sample(lambda.seq, 100))
    } else{
      lambdas <- sort(lambda.seq)
    }
  }
  
  names(lambdas) <- lambdas 

  # Shuffle data
  if(!use.bootstrap){
    if(sort) data <- data[sample(1:nrow(data), size = nrow(data)),]
    folds <- cut(seq(1,nrow(data)), breaks = nfolds, labels = FALSE)
  } else{
    if(sort) data <- data[sample(1:nrow(data), size = nrow(data)),]
    folds <- lapply(1:nfolds, function(i) unique(sample(1:nrow(data), nrow(data), replace = TRUE)))
  }
  
  in.train <- in.test <- trees <- list()
  
  # divide samples into training and testing based on n.folds specified
  for(k in 1:nfolds){
    if(!use.bootstrap){ # use traditional CV
      in.train[[k]] <- data[-which(folds==k,arr.ind=TRUE),]
      in.test[[k]]  <- data[which(folds==k,arr.ind=TRUE),]
    } else{
      in.train[[k]] <- data[folds[[k]],]
      in.test[[k]] <- data[-folds[[k]],]
    }
  }
  
  pboptions(type = ifelse(verbose, "timer", "none"))
  out.lambdas <- pblapply(lambdas, function(lam){

    # Grow initial tree to get sequence of alpha values
    tre <- rcDT(data = data, 
                split.var = split.var, 
                test = NULL, 
                min.ndsz = N0, 
                n0 = n0, 
                efficacy = efficacy, 
                risk = risk, 
                col.trt = col.trt,
                col.prtx = col.prtx,
                lambda = lam,
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
    
    full.tre.prune <- prune(tre, 
                            a = 0, 
                            ctgs = ctgs,
                            risk.control = risk.control, 
                            risk.threshold = risk.threshold, 
                            lambda = lam)
    
    alphas <- c(0, do.call(c, lapply(1:(nrow(full.tre.prune$result)-1), function(rtt){
      out.alpha <- (as.numeric(full.tre.prune$result$V[rtt]) - as.numeric(full.tre.prune$result$V[rtt+1])) / 
        (as.numeric(full.tre.prune$result$size.tmnl[rtt]) - as.numeric(full.tre.prune$result$size.tmnl[rtt+1]))
      return(out.alpha)})))
    names(alphas) <- alphas
    
    # Grow cross validation sample trees
    trees <- lapply(1:nfolds, function(n) 
      rcDT(data = in.train[[n]], 
           test = in.test[[n]], 
           split.var = split.var, 
           risk.control = risk.control, 
           risk.threshold = risk.threshold, 
           lambda = lam,
           efficacy = efficacy,
           risk = risk, 
           col.trt = col.trt,
           col.prtx = col.prtx,
           min.ndsz = N0, 
           n0 = n0, 
           AIPWE = AIPWE,
           ctg = ctg, 
           mtry = mtry,
           stabilize = stabilize, 
           stabilize.type = stabilize.type,
           max.depth = max.depth, 
           use.other.nodes = use.other.nodes))

    pruned.trees.alpha <- lapply(1:nfolds, function(tt){
      if(!is.null(dim(trees[[tt]]$tree) & nrow(trees[[tt]]$tree) != 1)){
        tmp <- prune(trees[[tt]], 
                     a = 0, 
                     test = trees[[tt]]$test, 
                     AIPWE = FALSE, 
                     ctgs = ctgs, 
                     risk.control = risk.control, 
                     risk.threshold = risk.threshold, 
                     lambda = lam)
        pruned.alphas <- do.call(c, lapply(alphas, function(aaa){
          idx <- which.max(as.numeric(tmp$result$V) - aaa * (as.numeric(tmp$result$size.tmnl)))
          return(as.numeric(tmp$result$V.test[idx]))
        }))
        return(list(V.test = pruned.alphas, 
                    pruned = tmp))
      }
    })
    
    best.alpha <- alphas[which.max(apply(do.call(rbind, lapply(pruned.trees.alpha, "[[", "V.test")), 2, mean, na.rm = TRUE))]
    best.tree.idx <- which.max(as.numeric(full.tre.prune$result$V) - best.alpha * as.numeric(full.tre.prune$result$size.tmnl))
    best.tree <- full.tre.prune$subtrees[[best.tree.idx]]
    
    best.tree.summary <- apply(do.call(rbind, lapply(1:nfolds, function(n){
      pruned <- pruned.trees.alpha[[n]]$pruned
      best.subtre.idx <- which.max(as.numeric(pruned$result$V) - best.alpha * as.numeric(pruned$result$size.tmnl))
      if(best.subtre.idx > length(pruned$subtrees)){
        return(c(Train.Benefit = sum(trees[[n]]$y / in.train[[n]][,col.prtx]) / 
                   sum(1 / in.train[[n]][,col.prtx]),
                 Train.Risk = sum(in.train[[n]][,risk] / in.train[[n]][,col.prtx]) / 
                   sum(1 / in.train[[n]][,col.prtx]),
                 Test.Benefit = sum(trees[[n]]$y.test / in.test[[n]][,col.prtx]) / 
                   sum(1 / in.test[[n]][,col.prtx]),
                 Test.Risk = sum(in.test[[n]][,risk] / in.test[[n]][,col.prtx]) / 
                   sum(1 / in.test[[n]][,col.prtx])))
      } else{
        tmp.tre <- pruned$subtrees[[best.subtre.idx]]
        
        pr.train <- predict.ITR(tmp.tre, in.train[[n]], split.var)$trt.pred
        pr.test <-  predict.ITR(tmp.tre, in.test[[n]], split.var)$trt.pred
        return(c(Train.Benefit = sum(trees[[n]]$y * (in.train[[n]][,col.trt] == pr.train) / in.train[[n]][,col.prtx]) / 
                   sum((in.train[[n]][,col.trt] == pr.train) / in.train[[n]][,col.prtx]),
                 Train.Risk = sum(in.train[[n]][,risk] * (in.train[[n]][,col.trt] == pr.train) / in.train[[n]][,col.prtx]) / 
                   sum((in.train[[n]][,col.trt] == pr.train) / in.train[[n]][,col.prtx]),
                 Test.Benefit = sum(trees[[n]]$y.test * (in.test[[n]][,col.trt] == pr.test) / in.test[[n]][,col.prtx]) / 
                   sum((in.test[[n]][,col.trt] == pr.test) / in.test[[n]][,col.prtx]),
                 Test.Risk = sum(in.test[[n]][,risk] * (in.test[[n]][,col.trt] == pr.test) / in.test[[n]][,col.prtx]) / 
                   sum((in.test[[n]][,col.trt] == pr.test) / in.test[[n]][,col.prtx])))
      }
    })), 2, mean, na.rm = TRUE)
    
    return(list(summaries = c(alpha = best.alpha, 
                              lambda = lam),
                best.tree = best.tree,
                full.tree = tre,
                full.tree.pruned = full.tre.prune,
                best.tree.summary = best.tree.summary))
  })


  tree.summaries <- do.call(rbind, lapply(out.lambdas, "[[", "best.tree.summary"))
  
  idx <- which(tree.summaries[,"Test.Risk"] < 1.01*risk.threshold)
  idx2 <- which.max(tree.summaries[idx,"Test.Benefit"])
  best.tree <- out.lambdas[idx][[idx2]]$best.tree
  best.lambda <- lambdas[idx][idx2]
  best.alpha <- out.lambdas[idx][[idx2]]$summaries[1]
  
  end.time <- Sys.time()
  
  setNames(list(best.tree,
                best.alpha,
                best.lambda,
                out.lambdas[idx][[idx2]]$full.tree, 
                out.lambdas[idx][[idx2]]$full.tree.pruned$result, 
                out.lambdas[idx][[idx2]]$full.tree.pruned$subtrees,
                tree.summaries,
                in.train, 
                in.test,
                end.time - start.time), 
           c("best.tree",
             "alpha", 
             "lambda", 
             "full.tree", 
             "pruned.tree",
             "subtrees", 
             "best.tree.summaries",
             "in.train", 
             "in.test",
             "elapsed.time"))
}

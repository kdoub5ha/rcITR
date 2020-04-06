\name{rcITR-package}
\alias{rcITR-package}
\alias{rcITR}
\docType{package}
\title{
\packageTitle{rcITR}
}
\description{
\packageDescription{rcITR}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{rcITR}
\packageIndices{rcITR}
  Tool for discovery of risk controlled individualized treatment rules
}
\author{
\packageAuthor{rcITR}

Maintainer: \packageMaintainer{rcITR}
}
\references{
[In Submission] Doubleday, K., Zhou, J., Fu, H. (2020), "Risk Controlled Decision Trees and Random Forests for Individualized Treatment Recommendations".
}
\keyword{ package }

\examples{
  \dontrun{
    set.seed(123)
    # Generate example data
    dat <- generateData(n = 1000)
    
    # Construct rcDT
    tre1 <- rcDT(data = dat, 
                 split.var = 1:10,
                 risk.threshold = 2.75,
                 lambda = 1,
                 efficacy = "y",
                 risk = "r",
                 col.trt = "trt",
                 col.prtx = "prtx")
    
    # Prune and select rcDT
    pruned1 <- prune(tre1, a = 0, risk.threshold = 2.75, lambda = 1)
    rcDT.fit <- rcDT.select(dat = dat, 
                            split.var = 1:10,
                            lambda = 1,
                            risk.threshold = 2.75, 
                            efficacy = "y",
                            risk = "r",
                            col.trt = "trt",
                            col.prtx = "prtx",
                            nfolds = 5)
    
    # Construct rcRF
    set.seed(2)
    rcRF.fit <- rcRF(dat, 
                     split.var = 1:10, 
                     efficacy = "y", 
                     risk = "r",
                     col.trt = "trt",
                     col.prtx = "prtx",
                     risk.threshold = 2.75,
                     ntree = 100, 
                     lambda = 0.5)
    
    # Variable importances
    VI <- Variable.Importance.ITR(rcRF.fit, sort = FALSE)
    
    # Predictions
    preds.rcDT <- predict.ITR(rcDT.fit$best.tree.alpha, 
                              new.data = dat, 
                              split.var = 1:10)
    preds.rcRF <- predict.ITR(rcRF.fit, 
                              new.data = dat, 
                              split.var = 1:10)
  }
}
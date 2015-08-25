library(partykit)
library(lme4)
library(Formula)
load("datasets1")
load("descriptions1")

lmertree <- function(lmtreeformula, randomformula, data, #subset = NULL, 
                       initialRandomEffects = NULL, 
                       ErrorTolerance = 0.001, MaxIterations = 1000, 
                       verbose = T, plotting=T,
                       lmer.control = lmerControl(), ...)
{
  iteration <- 0
  if (!is.null(initialRandomEffects)) {data$offset <- initialRandomEffects}
  if (is.null(initialRandomEffects)) {data$offset <- rep(0, times=dim(data)[1])}
  continuecondition <- T
  oldloglik <- -Inf
  while (continuecondition) {
    iteration <- iteration+1
    tree <- lmtree(lmtreeformula, data, offset=offset, ...)
    if(verbose){plot(tree)}
    if(length(tree)>1) {
      data$treeresponse <- predict(tree, newdata=data)
    }
    if(length(tree)==1) { # if treedepth=1, predict.lmtree doesn't work, then use unpartitioned lm
      lmformula <- formula(as.Formula(lmtreeformula), lhs=1, rhs = -2)
      data$treeresponse <- predict(lm(lmformula, data=data), newdata=data)
    }
    lme <- lmer(randomformula, data=data, offset=treeresponse)
    data$offset <- predict(lme, newdata=data)
    newloglik <- logLik(lme)    
    continuecondition <- (newloglik - oldloglik > ErrorTolerance) & (iteration < MaxIterations) 
    oldloglik <- newloglik
    if(verbose){print(newloglik)}
  } 
  result <- list(Tree = tree, RandomEffectModel = lme, RandomEffects = ranef(lme), 
                 BetweenMatrix = VarCorr(lme), ErrorVariance = attr(VarCorr(lme),"sc")^2, 
                 data = data, logLik = newloglik, IterationsUsed = iteration, 
                 MaxIterations=MaxIterations, initialRandomEffects = initialRandomEffects, 
                 TreeFormula = lmtreeformula, RandomFormula = randomformula, #Subset = subs, 
                 ErrorTolerance = ErrorTolerance, #residuals = residuals, 
                 tree.control = list(...), lmer.control = lmer.control)
}



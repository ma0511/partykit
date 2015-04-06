lmertree <- function(lmtreeformula, randomformula, data, #subset = NULL, 
  initialRandomEffects = NULL, 
  ErrorTolerance = 0.001, MaxIterations = 1000, 
  verbose = TRUE, plotting = FALSE,
  lmer.control = lmerControl(), ...)
{
  iteration <- 0L
  if (!is.null(initialRandomEffects)) { data$offset <- initialRandomEffects }
  if (is.null(initialRandomEffects)) { data$offset <- rep(0, times=dim(data)[1L]) }
  continuecondition <- TRUE
  oldloglik <- -Inf
  while (continuecondition) {
    iteration <- iteration + 1L
    tree <- lmtree(lmtreeformula, data, offset = offset, ...)
    if(plotting){plot(tree)}
    if(length(tree) > 1L) {
      data$treeresponse <- treeresponse <- predict(tree, newdata = data)
    }
    if(length(tree) == 1L) { # if treedepth=1, predict.lmtree doesn't work, then use unpartitioned lm
      lmformula <- formula(as.Formula(lmtreeformula), lhs = 1L, rhs = -2L)
      data$treeresponse <- treeresponse <- predict(lm(lmformula, data = data), newdata = data)
    }
    lme <- lmer(randomformula, data = data, offset = treeresponse)
    data$offset <- predict(lme, newdata = data)
    newloglik <- logLik(lme)    
    continuecondition <- (newloglik - oldloglik > ErrorTolerance) & (iteration < MaxIterations) 
    oldloglik <- newloglik
    if(verbose) {
      print(newloglik)
    }
  } 
  result <- list(
    Tree = tree,
    RandomEffectModel = lme,
    RandomEffects = ranef(lme), 
    BetweenMatrix = VarCorr(lme),
    ErrorVariance = attr(VarCorr(lme),"sc")^2, 
    data = data,
    logLik = newloglik,
    IterationsUsed = iteration, 
    MaxIterations=MaxIterations,
    initialRandomEffects = initialRandomEffects, 
    TreeFormula = lmtreeformula,
    RandomFormula = randomformula,
    #Subset = subs, 
    ErrorTolerance = ErrorTolerance,
    #residuals = residuals, 
    tree.control = list(...),
    lmer.control = lmer.control
  )
}

glmertree <- function(glmtreeformula, randomformula, data, #subset = NULL, 
  family = "binomial",
  initialRandomEffects = NULL, 
  ErrorTolerance = 0.001, MaxIterations = 1000, 
  verbose = TRUE, plotting = FALSE,
  glmer.control = glmerControl(), ...)
{
  iteration <- 0L
  if (!is.null(initialRandomEffects)) { data$offset <- initialRandomEffects }
  if (is.null(initialRandomEffects)) { data$offset <- rep(0, times = dim(data)[1L]) }
  continuecondition <- TRUE
  oldloglik <- -Inf
  while (continuecondition) {
    iteration <- iteration + 1L
    tree <- glmtree(glmtreeformula, data, family = family, offset = offset, ...)
    if(plotting){plot(tree)}
    if(length(tree) > 1L) {
      data$treeresponse <- treeresponse <- predict(tree, newdata = data, type = "link")
    }
    if(length(tree) == 1L) { # if treedepth=1, predict.glmtree doesn't work, then use unpartitioned glm
      glmformula <- formula(as.Formula(glmtreeformula), lhs = 1L, rhs = -2L)
      data$treeresponse <- treeresponse <- predict(glm(glmformula, family = family, data = data), newdata = data, type = "link")
    }
    glme <- glmer(randomformula, family = family, data = data, offset = treeresponse)
    data$offset <- predict(glme, newdata = data, type = "link")
    newloglik <- logLik(glme)    
    continuecondition <- (newloglik - oldloglik > ErrorTolerance) & (iteration < MaxIterations) 
    oldloglik <- newloglik
    if(verbose) {
      print(newloglik)
    }
  } 
  result <- list(
    Tree = tree,
    RandomEffectModel = glme,
    RandomEffects = ranef(glme), 
    BetweenMatrix = VarCorr(glme),
    ErrorVariance = attr(VarCorr(glme),"sc")^2, 
    data = data,
    logLik = newloglik,
    IterationsUsed = iteration, 
    MaxIterations = MaxIterations,
    initialRandomEffects = initialRandomEffects, 
    TreeFormula = glmtreeformula,
    RandomFormula = randomformula,
    #Subset = subs, 
    ErrorTolerance = ErrorTolerance,
    #residuals = data$residuals, 
    tree.control = list(...),
    glmer.control = glmer.control
  )
}

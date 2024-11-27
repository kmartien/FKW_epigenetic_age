
# standard 10-fold cv.glmnet run at given alpha
CVglmnet10K <- function(df, sites, age.transform, alpha) {
  if (age.transform == 'ln') df$age.best <- log(df$age.best)
  cv.glmnet(
    x = as.matrix(df[, sites]),
    y = df$age.best,
    weights = df$wt,
    alpha = alpha,
  ) 
}

# return glmnet median MAE at minimum lambda for alpha and nrep replicates
median.MAE.min <- function(alpha, df, sites, age.transform, nrep) {
  parallel::mclapply(1:nrep, function(i) {
    cv.fit <- tryCatch({
      CVglmnet10K(df, sites, age.transform, alpha)
    }, error = function(e) NULL)
    if(is.null(cv.fit)) NA else {
      predicted.age <- predict(cv.fit, as.matrix(df[,sites]), s = 'lambda.min')
      if(age.transform == 'ln') predicted.age <- exp(predicted.age)
      median(abs(df$age.best - exp(predict(cv.fit, as.matrix(df[,sites]), s = 'lambda.min'))))
    } 
  }, mc.cores = 6) |> 
    unlist() |> 
    median(na.rm = TRUE)
}

#' run randomForest over sampsize and mtry grid and report deviance stats
rf.param.grid.search <- function(df, sites, age.transform, ntree){
  # grid of sampsize and mtry to optimize MSE over
  params <- expand.grid(
    sampsize = seq(round(nrow(df)*.3), round(nrow(df)*.7), by = 1), 
    mtry = seq(round(length(sites)*.2), round(length(sites)*.5), by = 1), 
    KEEP.OUT.ATTRS = FALSE
  )
  
  parallel::mclapply(1:nrow(params), function(i){
    if(age.transform == 'ln') df$age.best <- log(df$age.best)
    rf <- randomForest(
      y = df$age.best,
      x = df[, sites.to.keep],
      mtry = params$mtry[i],
      ntree = ntree,
      sampsize = params$sampsize[i],
      replace = FALSE
    )
    data.frame( 
      sampsize = params$sampsize[i],
      mtry = params$mtry[i],
      mse = rf$mse[length(rf$mse)],
      rsq = rf$rsq[length(rf$rsq)],
      pct.var = 100 * rf$rsq[length(rf$rsq)]
    )
  }, mc.cores = 6) |> 
    bind_rows() |> 
    as.data.frame()
}


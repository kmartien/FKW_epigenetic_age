unlink('~/temp', recursive = TRUE, force = TRUE)
dir.create('~/temp')
Sys.setenv(TEMP = '~/temp', TMP = '~/temp')

library(tidyverse)

# estimate skew-normal parameters from best, min, and max
sn.params <- function(
    age.mode, age.min, age.max, p, shape.const = 5, scale.const = 0.5, 
    maxit = 5000
) {
  age.mean <- mean(c(age.min, age.max))
  shape <- shape.const * (age.mean - age.mode) / (age.mean - age.min)
  scale <- scale.const * (age.max - age.min)
  optim(
    par = c(age.mode, scale, shape), 
    fn = function(dp, obs, p) {
      if(dp[2] <= 0) dp[2] <- .Machine$double.eps 
      abs(obs[1] - swfscMisc::sn.mode(dp)) +
        abs(obs[2] - sn::qsn(1 - p, dp = dp, solver = 'RFB')) +
        abs(obs[3] - sn::qsn(p, dp = dp, solver = 'RFB'))
    }, 
    obs = c(age.mode, age.min, age.max), 
    p = p,
    control = list(
      maxit = maxit, 
      abstol = .Machine$double.eps,
      reltol = 1e-12
    )
  )
}

cwsnDensity <- function(age, age.range, dp, wt) {
  dens.sn <- sn::dsn(age, dp = dp)
  wtd.unif <- (1 - wt) / age.range
  (wtd.unif + (dens.sn * wt))
}

maxCWSNdensity <- function(df) {
  max(sapply(1:nrow(df), function(i) {
    x <- df[i, ]
    cwsnDensity(
      x$age.best, 
      x$age.range,
      dp = c(x$sn.location, x$sn.scale, x$sn.shape),
      wt = x$confidence.wt
    )
  }))
}

rcwsn <- function(n, age.min, age.max, dp, wt) {
  accepted <- numeric(0)
  max.lik <- cwsnDensity(swfscMisc::sn.mode(dp), age.max - age.min, dp, wt)
  while(length(accepted) < n) {
    unif.age <- runif(1, age.min, age.max) 
    if(runif(1, 0, 1) < cwsnDensity(unif.age, age.max - age.min, dp, wt) / max.lik) {
      accepted <- c(accepted, unif.age)
    }
  }
  accepted
}

ranAges <- function(df) {
  sapply(1:nrow(df), function(i) {
    rcwsn(
      n = 1, 
      age.min = df$age.min[i], 
      age.max = df$age.max[i], 
      dp = unlist(df[i, c('sn.location', 'sn.scale', 'sn.shape')]),
      wt = df$confidence.wt[i]
    )
  })
}

crcAgeDist <- function(
    id, age.df,
    n.samples = if(match.arg(type) == 'samples') 10000 else 1000,
    type = c('samples', 'density'),
    from = 0, to = 80
) {
  type <- match.arg(type)
  x <- age.df |> 
    filter(swfsc.id == id) |> 
    as.list()
  
  age <- if(type == 'samples') {
    runif(n.samples * 10, x$age.min, x$age.max)
  } else {
    seq(from, to, length.out = n.samples)
  }
  
  cwsn.dens <- cwsnDensity(
    age = age, 
    age.range = x$age.range, 
    dp = c(x$sn.location, x$sn.scale, x$sn.shape), 
    wt = x$confidence.wt
  )
  
  if(type == 'samples') {
    sample(age[cwsn.dens > 0], n.samples, p = cwsn.dens[cwsn.dens > 0])
  } else {
    data.frame(swfsc.id = id, age = age, density = cwsn.dens) 
  }
}

ranMeth <- function(cpg) {
  meth <- cpg |> 
    mutate(uncorrected.ran.meth = (rbinom(n(), corrected.cov, pct.meth) / corrected.cov),
           corrected.ran.meth = 1 - ((1 - uncorrected.ran.meth) / conversion)) |> 
    select(id, loc.site, corrected.ran.meth)
  loc.min.non.zero <- filter(meth, corrected.ran.meth > 0) |> 
    group_by(loc.site) |> 
    summarise(min = min(corrected.ran.meth, na.rm = TRUE))
  left_join(meth, loc.min.non.zero, by = 'loc.site') |> 
    mutate(ran.meth = ifelse(corrected.ran.meth <= 0, (min / 2), corrected.ran.meth),
           ran.meth = swfscMisc::logOdds(ran.meth))
}

sampleAgeMeth <- function(ages, meth) {
  ages |>  
    mutate(age.ran = ranAges(ages)) |> 
    left_join(
      ranMeth(meth) |> 
        rename(swfsc.id = id) |> 
        select(swfsc.id, loc.site, ran.meth) |>
        pivot_wider(names_from = 'loc.site', values_from = 'ran.meth'),
      by = 'swfsc.id'
    ) 
}

calc.ci.wt <- function(model.df, ndraws = 1000){
  left_join(
    model.df |> 
      mutate(age.confidence = factor(age.confidence)), 
    
    bind_rows(
      lapply(1:ndraws, function(i){
        mutate(model.df, age.ran = ranAges(model.df)) |> 
          select(swfsc.id, age.ran)
      })) |> 
      group_by(swfsc.id) |> 
      summarise(lci = quantile(age.ran, probs = c(0.25)),
                uci = quantile(age.ran, probs = c(0.75)),
                ci.range = uci - lci) |> 
      select(swfsc.id, ci.range)
  ) |> 
    mutate(ci.wt = 1 - (ci.range / age.range),
           ci.wt = (ci.wt - 0.3) / 0.4) |> 
    pull(ci.wt)
}

selectCpGsites <- function(sites.2.use, minCR, weight){
  if(sites.2.use == 'RFsites'){
    # select important sites from Random Forest tuned with all samples
    chosen.sites <- readRDS('R/rf_tuning/rf_chosen_sites.rds') |>
      filter(cr == minCR & wt == weight & pval <= 0.05) |>
      pull('loc.site')
  }
  if(sites.2.use == 'glmnet'){
    # select chosen sites from glmnet tuned with all samples, alpha = 0.5
    chosen.sites <- readRDS('R/glmnet_tuning/glmnet.chosen.sites.rds') |> 
      filter(cr == minCR & wt == weight & count >= 500 & sites != '(Intercept)') |> 
      pull('sites')
  }
  return(chosen.sites)
}

# fit GAM model to training data
fitTrainGAM <- function(train.df, sites, resp, age.transform) {
  if (age.transform == 'ln') train.df[[resp]] <- log(train.df[[resp]])
  # run GAM starting at k = 3 until convergence
  fit <- NULL
  for(k in 3:10) {
    fit <- tryCatch(
      gam(
        formula = as.formula(paste(
          resp, '~',
          paste0('te(', sites, ', k = ', k, ', bs = "ts")', collapse = ' + ')
        )),
        data = train.df,
        weights = train.df$wt,
        select = TRUE
      ),
      error = function(e) NULL
    )
    if(!is.null(fit)) break
  }
  if(is.null(fit)) NULL else list(gam = fit, k = k)
}

# predict testing data from fitted model
predictTestGAM <- function(fit, test.df, resp, age.transform){
  if(is.null(fit)) return(NULL)
  
  pred <- predict(
    fit$gam,
    newdata = test.df, 
    type = 'response',
    se.fit = TRUE
  )

  tibble(
    swfsc.id = test.df$swfsc.id,
    k = fit$k,
    age.resp = test.df[[resp]],
    age.pred = unname(pred$fit),
    lci = unname(pred$fit + (pred$se.fit * qnorm(0.025))),
    uci = unname(pred$fit + (pred$se.fit * qnorm(0.975)))
  ) |> 
    mutate(
      age.pred = if (age.transform == 'ln') exp(age.pred) else age.pred,
      age.pred = ifelse(age.pred < 0, 0, age.pred),
      age.pred = ifelse(age.pred > 80, 80, age.pred),
      lci = if (age.transform == 'ln') exp(lci) else lci,
      lci = ifelse(lci < 0, 0, lci),
      lci = ifelse(lci > 80, 80, lci),
      uci = if (age.transform == 'ln') exp(uci) else uci,
      uci = ifelse(uci < 0, 0, uci),
      uci = ifelse(uci > 80, 80, uci),
      resid = age.pred - age.resp,
      dev = abs(resid),
      sq.err = resid ^ 2
    )
}

predictAllIDsGAM <- function(train.df, model.df, sites, resp, age.transform = 'none') {
  rbind( 
    # cross-validation model for training samples
    lapply(train.df$swfsc.id, function(cv.id) {
      fitTrainGAM(filter(train.df, swfsc.id != cv.id), sites, resp, age.transform) |> 
        predictTestGAM(filter(model.df, swfsc.id == cv.id), resp, age.transform)
    }) |> 
      bind_rows(),
    # full model to predict test samples
    fitTrainGAM(train.df, sites, resp, age.transform) |> 
      predictTestGAM(filter(model.df, !swfsc.id %in% train.df$swfsc.id), resp, age.transform)
  ) 
}

fitTrainSVM <- function(df, sites, resp, svm.params, age.transform) {
  if (age.transform == 'ln') df[[resp]] <- log(df[[resp]])
  fit <- svm(
    formula = as.formula(paste0(resp, ' ~ .')), 
    data = select(df, c(resp, all_of(sites))),
    cost = svm.params$cost,
    gamma = svm.params$gamma)
}


predictTestSVM <- function(fit, cv.df, sites, resp, age.transform){
  if(is.null(fit)) return(NULL)
  
  pred <- predict(
    fit, 
    select(cv.df, all_of(sites))
  )
  if (age.transform == 'ln') pred <- exp(pred)
  
  tibble(
    swfsc.id = cv.df$swfsc.id,
    age.resp = cv.df[[resp]],
    age.pred = unname(pred),
  ) |> 
    mutate(
      age.pred = ifelse(age.pred < 0, 0, age.pred),
      age.pred = ifelse(age.pred > 80, 80, age.pred),
      resid = age.pred - age.resp,
      dev = abs(resid),
      sq.err = resid ^ 2
    )
}

predictAllIDsSVM <- function(train.df, model.df, sites, resp, svm.params, age.transform = 'none') {
  rbind( 
    # cross-validation model for training samples
    lapply(train.df$swfsc.id, function(cv.id) {
      fitTrainSVM(filter(train.df, swfsc.id != cv.id), sites, resp, svm.params, age.transform) |> 
        predictTestSVM(filter(model.df, swfsc.id == cv.id), sites, resp, age.transform)
    }) |> 
      bind_rows(),
    # full model to predict test samples
    fitTrainSVM(train.df, sites, resp, svm.params,age.transform) |> 
      predictTestSVM(filter(model.df, !swfsc.id %in% train.df$swfsc.id), sites, resp, age.transform)
  ) 
}

fitTrainRF <- function(df, sites, resp, rf.params, age.transform) {
  if (age.transform == 'ln') df[[resp]] <- log(df[[resp]])
  fit <- randomForest(
    formula = as.formula(paste0(resp, ' ~ .')), 
    data = select(df, c(resp, all_of(sites))),
    weights = df$wt,
    mtry = rf.params$mtry,
    ntree = 10000,
    sampsize = rf.params$sampsize,
    replace = FALSE
  )
}

predictTestRF <- function(fit, test.df, sites, resp, age.transform){

  pred <- 
    if(is.null(fit)){
      fitTrainRF(test.df, sites, resp, rf.params, age.transform)$predicted
    } else {
      predict(
        fit, 
        select(test.df, all_of(sites))
      )
    }
  if (age.transform == 'ln') pred <- exp(pred)
      
  tibble(
    swfsc.id = test.df$swfsc.id,
    age.resp = test.df[[resp]],
    age.pred = unname(pred),
  ) |> 
    mutate(
      age.pred = ifelse(age.pred < 0, 0, age.pred),
      age.pred = ifelse(age.pred > 80, 80, age.pred),
      resid = age.pred - age.resp,
      dev = abs(resid),
      sq.err = resid ^ 2
    )
}

predictAllIDsRF <- function(train.df, model.df, sites, resp, rf.params, age.transform = 'none') {
  rbind( 
    # OOB predictions for training samples
        predictTestRF(fit = NULL, train.df, sites, resp, age.transform),
    # full model to predict test samples
    fitTrainRF(train.df, sites, resp, rf.params, age.transform) |> 
      predictTestSVM(filter(model.df, !swfsc.id %in% train.df$swfsc.id), sites, resp, age.transform)
  ) 
}
  
fitTrainENR <- function(df, sites, resp, alpha, age.transform) {
  if (age.transform == 'ln') df[[resp]] <- log(df[[resp]])
  fit <- cv.glmnet(
    x = as.matrix(df[, sites]),
    y = df[,resp],
    weights = df$wt,
    alpha = alpha
  ) 
}

predictTestENR <- function(fit, cv.df, sites, resp, age.transform){
  if(is.null(fit)) return(NULL)
  
  pred <- predict(
    fit, 
    as.matrix(cv.df[, sites]),
    s = 'lambda.min'
  )
  if (age.transform == 'ln') pred <- exp(pred)
    
  tibble(
    swfsc.id = cv.df$swfsc.id,
    age.resp = cv.df[[resp]],
    age.pred = pred[1],
  ) |> 
    mutate(
      age.pred = ifelse(age.pred < 0, 0, age.pred),
      age.pred = ifelse(age.pred > 80, 80, age.pred),
      resid = age.pred - age.resp,
      dev = abs(resid),
      sq.err = resid ^ 2
    )
}

predictAllIDsENR <- function(train.df, model.df, sites, resp, alpha, age.transform = 'none') {
  rbind( 
    # cross-validation model for training samples
    lapply(train.df$swfsc.id, function(cv.id) {
      fitTrainENR(filter(train.df, swfsc.id != cv.id), sites, resp, alpha, age.transform) |> 
        predictTestENR(filter(model.df, swfsc.id == cv.id), sites, resp, age.transform)
    }) |> 
      bind_rows(),
    # full model to predict test samples
    fitTrainENR(train.df, sites, resp, alpha, age.transform) |> 
      predictTestENR(filter(model.df, !swfsc.id %in% train.df$swfsc.id), sites, resp, age.transform)
  ) 
}



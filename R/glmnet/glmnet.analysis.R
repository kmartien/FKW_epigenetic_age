library(glmnet)
library(dplyr)
library(tidyverse)
source('R/misc_funcs.R')
source('R/tuning_funcs.R')
load("data/age_and_methylation_data.rdata")
load("data/model.params.rda")

nrep <- 1000
ncores <- 10

top.models <- readRDS('data/top.models.rds')
model.params <- filter(model.params, model %in% top.models$glmnet)

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

age.df$ci.wt <- calc.ci.wt(age.df)

lapply(1:nrow(model.params), function(i){
  print(i)
  print(date())
  params <- model.params[i,]
  description <- paste0("model", params$model)
  weight <- params$weight
  minCR <- params$training.cr
  sites.2.use <- params$site.select.regr.meth
  site.select.cr <- params$site.select.cr
  age.transform <-params$age.transform
  
  sites <- sites.to.keep
  if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use, site.select.cr, weight, age.transform)
  
  age.df <- age.df |> 
    mutate(
      wt = if(weight == 'ci.wt') ci.wt else {
        if (weight == 'CR') age.confidence else {
          if (weight == 'sn.wt') confidence.wt else 1
        }
      })
  
  model.df <- left_join(
    age.df,
    rownames_to_column(logit.meth, var = 'swfsc.id'),
    by = 'swfsc.id'
  )
  
  # age prediction
  train.df <- filter(model.df, age.confidence >= minCR)
  
  # find optimal alpha (lowest median cvm at minimum lambda)
  if(file.exists(file = paste0('R/glmnet_tuning/optim_alpha_', description, '.rds'))){
    opt.alpha <- readRDS(paste0('R/glmnet_tuning/optim_alpha_', description, '.rds'))
  } else {
    opt.alpha <- optim(
      par = c(alpha = 0.3),
      fn = median.cvm.min,
      method = "Brent",
      lower = 0.001,
      upper = 0.999,
      df = train.df,
      sites = sites,
      age.transform = age.transform,
      nrep = nrep
    )
    saveRDS(opt.alpha, file = paste0('R/glmnet_tuning/optim_alpha_', description, '.rds'))
  }
  
  # Best age and methylation estimates --------------------------------------
  
  # pred.best <- predictAllIDsENR(train.df, model.df, sites, 'age.best', opt.alpha$par, age.transform) #|>
  # saveRDS(pred.best, paste0('R/glmnet/glmnet_best_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))

  # Random age and random methylation estimates -----------------------------

  pred.ranAgeMeth <- parallel::mclapply(1:nrep, function(j) {
    # random sample of ages and methylation
    ran.df <- sampleAgeMeth(age.df, cpg)

    train.df <- filter(ran.df, age.confidence >= minCR)
    predictAllIDsENR(train.df, ran.df, sites, 'age.ran', opt.alpha$par, age.transform)
  }, mc.cores = ncores) |>
    bind_rows()
  saveRDS(pred.ranAgeMeth, paste0('R/glmnet/glmnet_ranAgeMeth_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))
})


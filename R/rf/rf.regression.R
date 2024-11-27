library(randomForest)
library(rfPermute)
library(dplyr)
library(tidyverse)
source('R/misc_funcs.R')
source('R/tuning_funcs.R')
load("data/age_and_methylation_data.rdata")
load("data/model.params.rda")

nrep <- 1000
ntree <- 5000
ncores <- 10

model.params <- filter(model.params, num.sites > 100)

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

age.df$ci.wt <- calc.ci.wt(age.df)

lapply(10:nrow(model.params), function(i){
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
  
  model.df <- age.df |> 
    mutate(
      wt = if(weight == 'ci.wt') ci.wt else {
        if (weight == 'CR') age.confidence else {
          if (weight == 'sn.wt') confidence.wt else 1
        }
      }) |>
    left_join(
      rownames_to_column(logit.meth, var = 'swfsc.id'),
      by = 'swfsc.id'
    )
  
  # age prediction
  train.df <- filter(model.df, age.confidence >= minCR)

  if(file.exists(file = paste0('R/rf_tuning/rf_optim_params_', description, '.rda'))){
    load(paste0('R/rf_tuning/rf_optim_params_', description, '.rda'))
  } else {
    param.df  <- rf.param.grid.search(train.df, sites, age.transform, ntree)
    save(param.df, file = paste0('R/rf_tuning/rf_optim_params_', description, '.grid.search.rda'))
    rf.params <- filter(param.df, mse == min(param.df$mse)) |>
      select(c(mtry, sampsize))
    save(rf.params, file = paste0('R/rf_tuning/rf_optim_params_', description, '.rda'))
  }
  
  # Best age and methylation estimates --------------------------------------
  
  print('Best')
  pred <- if(minCR == 2) {
    # OOB predictions for training samples
    predictTestRF(fit = NULL, train.df, sites, 'age.best', rf.params, age.transform)
  } else {
    predictAllIDsRF(train.df, model.df, sites, 'age.best', rf.params, age.transform)  
  }
  saveRDS(pred, paste0('R/rf/rf_best_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))
  
  # Random age and random methylation estimates -----------------------------

  # print('RanAgeMeth')
  # pred.ranAgeMeth <- parallel::mclapply(1:nrep, function(j) {
  #   # random sample of ages and methylation
  #   ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params)
  #   
  #   if(minCR == 2) {
  #     # OOB predictions for training samples
  #     predictTestRF(fit = NULL, ran.df, sites, 'age.ran', age.transform)
  #   } else {
  #     predictAllIDsRF(filter(ran.df, age.confidence >= minCR), ran.df, sites, 'age.ran', rf.params, age.transform)
  #   }
  # }, mc.cores = ncores) |>
  #   bind_rows()
  # saveRDS(pred.ranAgeMeth, paste0('R/rf/rf_ranAgeMeth_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))
})

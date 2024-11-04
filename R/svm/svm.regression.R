library(dplyr)
library(tidyverse)
library(e1071)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")
load("data/model.params.rda")

model.params <- model.params[-c(31:36),]

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

lapply(64:nrow(model.params), function(i){
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
  if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use, site.select.cr, weight)
  
  model.df <- age.df |> 
    mutate(
      wt = if(weight == 'inv.var') 1/age.var else {
        if (weight == 'CR') age.confidence else 1
      }) |> 
    left_join(
      rownames_to_column(logit.meth, var = 'swfsc.id'),
      by = 'swfsc.id'
    )
  
  # age prediction
  train.df <- filter(model.df, age.confidence >= minCR)
  
  print(date())
  tune.obj <- tune(svm, 
                   age.best ~ .,
                   data = select(train.df, c(age.best, all_of(sites))),
                   ranges = list(
                     cost = 10^(seq(-4, 5, 0.1)),
                     gamma = 10^(seq(-5, 4, 0.1))),
                   tunecontrol = tune.control(sampling = "cross"),
                   cross = 10)
  print(date())
  
  pred <- predictAllIDsSVM(train.df, model.df, sites, 'age.best', tune.obj$best.parameters, age.transform)
  saveRDS(pred, paste0('R/svm/svm_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))

  save(tune.obj, file = paste0("R/svm/svm_tuning_", description, ".rda"))

})

         
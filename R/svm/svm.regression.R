library(dplyr)
library(tidyverse)
library(e1071)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")
load("data/model.params.rda")

nrep <- 1000
ncores <- 10

top.models <- readRDS('data/top.models.rds')
model.params <- filter(model.params, model %in% top.models$svm &
                         model != 21 & model != 24)

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

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
  
  if(file.exists(file = paste0("R/svm/svm_tuning_", description, ".rda"))){
    load(paste0("R/svm/svm_tuning_", description, ".rda"))
  } else {
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
    
    save(tune.obj, file = paste0("R/svm/svm_tuning_", description, ".rda"))
  }
  # pred <- predictAllIDsSVM(train.df, model.df, sites, 'age.best', tune.obj$best.parameters, age.transform)
  # saveRDS(pred, paste0('R/svm/svm_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))
  # 
  # Random age and random methylation estimates -----------------------------
   pred.ranAgeMeth <- parallel::mclapply(1:nrep, function(j) {
    # random sample of ages and methylation
    ran.df <- sampleAgeMeth(age.df, cpg)

    train.df <- filter(ran.df, age.confidence >= minCR)
    predictAllIDsSVM(train.df, ran.df, sites, 'age.ran', tune.obj$best.parameters, age.transform)
  }, mc.cores = ncores) |>
    bind_rows()
  saveRDS(pred.ranAgeMeth, paste0('R/svm/svm_ranAgeMeth_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))
  

})

         
rm(list = ls())
library(tidyverse)
library(mgcv)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")
load("data/model.params.rda")

ncores <- 10

model.params <- filter(model.params, num.sites < 40 & age.transform == 'ln')

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
  if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use, site.select.cr, weight)
  
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
  
  pred <- if (minCR == 2){
    # LOO cross validation for all samples
    lapply(model.df$swfsc.id, function(cv.id) {
      fitTrainGAM(filter(model.df, swfsc.id != cv.id), sites, 'age.best', age.transform) |>
        predictTestGAM(filter(model.df, swfsc.id == cv.id), 'age.best', age.transform)
    })|> bind_rows()
  } else {predictAllIDsGAM(train.df, model.df, sites, 'age.best', age.transform)}
  saveRDS(pred, paste0('R/gam/gam_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))
})


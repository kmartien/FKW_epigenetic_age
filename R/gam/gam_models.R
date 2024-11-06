rm(list = ls())
library(tidyverse)
library(mgcv)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")
load("data/model.params.rda")

nrep <- 1000
ncores <- 10

model.params <- filter(model.params, model == 72)

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

  # pred <- if (minCR == 2){
  #   # LOO cross validation for all samples
  #   lapply(model.df$swfsc.id, function(cv.id) {
  #     fitTrainGAM(filter(model.df, swfsc.id != cv.id), sites, 'age.best', age.transform) |>
  #       predictTestGAM(filter(model.df, swfsc.id == cv.id), 'age.best', age.transform)
  #   })|> bind_rows()
  # } else {predictAllIDsGAM(train.df, model.df, sites, 'age.best', age.transform)}
  # saveRDS(pred, paste0('R/gam/gam_best_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))

  # Random age and random methylation estimates -----------------------------

  message(format(Sys.time()), ' : RanAgeMeth - All')
  pred.ranAgeMeth <- parallel::mclapply(1:nrep, function(j) {
    print(j)
    print(date())
    ran.df <- sampleAgeMeth(age.df, cpg)
    if (minCR == 2){
      # LOO cross validation for all samples
      lapply(model.df$swfsc.id, function(cv.id) {
        fitTrainGAM(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', age.transform) |>
          predictTestGAM(filter(ran.df, swfsc.id == cv.id), 'age.ran', age.transform)
      })|> bind_rows()
    } else {predictAllIDsGAM(filter(ran.df, age.confidence >= minCR), ran.df, sites, 'age.ran', age.transform)}
  }, mc.cores = ncores) |> bind_rows()
  saveRDS(pred.ranAgeMeth, paste0('R/gam/gam_ranAgeMeth_minCR', minCR, '_', sites.2.use,'_cr', site.select.cr, '_', age.transform, '_', weight, '.rds'))
})


library(tidyverse)
library(randomForest)
library(rfPermute)
source('R/misc_funcs.R')
load('R/rf_tuning/rf_optim_params.rda')
load("data/age_and_methylation_data.rdata")

ntree <- 10000
nrep <- 1000

age.df$ci.wt <- calc.ci.wt(age.df)

age.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>
  left_join(
    rownames_to_column(logit.meth, var = 'swfsc.id')
  ) |> 
  column_to_rownames('swfsc.id') 

#------------------------------------------------------------------------
RFsites.ln <- do.call(rbind, lapply(10:nrow(rf.params), function(i){
  print(i)
  print(date())
  minCR <- rf.params$CR[i]
  weight <- rf.params$wt[i]
  mtry <- rf.params$mtry[i]
  sampsize <- rf.params$sampsize[i]
  age.transform <- rf.params$age.transform[i]
  model.df <- age.df |> 
    mutate(
      wt = if(weight == 'ci.wt') ci.wt else {
        if (weight == 'CR') age.confidence else {
          if (weight == 'sn.wt') confidence.wt else 1
        }
      },
      if(age.transform == 'ln') age.best = log(age.best))
  train.df <- filter(model.df, age.confidence >= minCR)
  sites <- sites.to.keep
  # All samples
  # run rfPermute at sampsize and mtry with minimum MSE
  rp <- rfPermute(
    y = model.df$age.best,
    x = model.df[, sites.to.keep],
    weights = model.df$wt,
    mtry = mtry,
    ntree = ntree,
    sampsize = sampsize,
    importance = TRUE,
    replace = FALSE,
    num.rep = nrep,
    num.cores = 10
  ) |> 
    importance() |> 
    as.data.frame() |> 
    rownames_to_column('loc.site') |> 
    rename(
      incMSE = '%IncMSE',
      pval = '%IncMSE.pval'
    ) |> 
    select(loc.site, incMSE, pval) 
  return(bind_cols(rp, cr = minCR, wt = weight))
}))

saveRDS(RFsites.ln, 'R/rf_tuning/rf_chosen_sites.ln.rds')

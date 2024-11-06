rm(list = ls())
library(tidyverse)
library(randomForest)
library(rfPermute)
source('R/misc_funcs.R')
source('R/tuning_funcs.R')
load("data/age_and_methylation_data.rdata")

ntree <- 1000
weight.schemes <- c('none', 'CR', 'ci.wt')
site.selection <- c('Allsites') #c('Allsites', 'RFsites', 'ENRsites')
age.transform <- c('ln')

age.df$ci.wt <- calc.ci.wt(age.df)

#' #' run randomForest over sampsize and mtry grid and report deviance stats
#' rf.param.grid.search <- function(df, sites, age.transform, ntree){
#'   # grid of sampsize and mtry to optimize MSE over
#'   params <- expand.grid(
#'     sampsize = seq(round(nrow(model.df)*.3), round(nrow(model.df)*.7), by = 1), 
#'     mtry = seq(round(length(sites)*.2), round(length(sites)*.5), by = 1), 
#'     KEEP.OUT.ATTRS = FALSE
#'   )
#'   
#'   parallel::mclapply(1:nrow(params), function(i){
#'     if(age.transform == 'ln') df$age.best <- log(df$age.best)
#'     rf <- randomForest(
#'       y = df$age.best,
#'       x = df[, sites.to.keep],
#'       mtry = params$mtry[i],
#'       ntree = ntree,
#'       sampsize = params$sampsize[i],
#'       replace = FALSE
#'     )
#'     data.frame( 
#'       sampsize = params$sampsize[i],
#'       mtry = params$mtry[i],
#'       mse = rf$mse[length(rf$mse)],
#'       rsq = rf$rsq[length(rf$rsq)],
#'       pct.var = 100 * rf$rsq[length(rf$rsq)]
#'     )
#'   }, mc.cores = 6) |> 
#'     bind_rows() |> 
#'     as.data.frame()
#' }
#' 
model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>
  left_join(
    rownames_to_column(logit.meth, var = 'swfsc.id')
  ) |> 
  column_to_rownames('swfsc.id')

rf.params.ln <- do.call(rbind, lapply(c(2, 3, 4), function(minCR){
  do.call(rbind, lapply(weight.schemes, function(weight){
    do.call(rbind, lapply(age.transform, function(age){
    print(paste0('minCR = ', minCR))
    print(paste0('weight = ', weight))
    print(paste0('age.transform = ', age))
    print(date())
    model.df <- model.df |> 
      mutate(
        wt = if(weight == 'ci.wt') ci.wt else {
          if (weight == 'CR') age.confidence else {
            if (weight == 'sn.wt') confidence.wt else 1
          }
        })
    train.df <- filter(model.df, age.confidence >= minCR)
    do.call(rbind, lapply(site.selection, function(sites.2.use){
      sites <- sites.to.keep
      if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)
      param.df  <- rf.param.grid.search(model.df, sites, age, ntree)
      save(param.df, file = paste0('R/rf_tuning/rf_minCR', minCR, '_', sites.2.use, '_.grid.search.rda'))
      return(tibble(s = sites.2.use, wt = weight, CR = minCR, age.transform = age,
                        filter(param.df, mse == min(param.df$mse)) |>
                          select(c(mtry, sampsize))
                        )
             )
    }))
    }))
  }))
}))
save(rf.params.ln, file = paste0('R/rf_tuning/rf_optim_params.ln.rda'))

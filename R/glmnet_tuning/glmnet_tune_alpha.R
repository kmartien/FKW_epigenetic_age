library(tidyverse)
library(glmnet)
source('R/misc_funcs.R')
source('R/tuning_funcs.R')
load("data/age_and_methylation_data.rdata")

nrep <- 1000
weight.schemes <- c('none', 'CR', 'ci.wt')
site.selection <- c('Allsites') #c('Allsites', 'RFsites', 'ENRsites')
age.transform <- 'ln'

age.df$ci.wt <- calc.ci.wt(age.df)

model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |> 
  left_join(
    rownames_to_column(logit.meth, var = 'swfsc.id')
  ) |> 
  column_to_rownames('swfsc.id') 
  

# # standard 10-fold cv.glmnet run at given alpha
# CVglmnet10K <- function(df, sites, age.transform, alpha) {
#   if (age.transform == 'ln') df$age.best <- log(df$age.best)
#   cv.glmnet(
#     x = as.matrix(df[, sites]),
#     y = df$age.best,
#     weights = df$wt,
#     alpha = alpha,
#   ) 
# }
# 
# return median cvm at minimum lambda for alpha and nrep replicates
# median.cvm.min <- function(alpha, df, sites, age.transform, nrep) {
#   parallel::mclapply(1:nrep, function(i) {
#     cv.fit <- tryCatch({
#       CVglmnet10K(df, sites, age.transform, alpha)
#     }, error = function(e) NULL)
#     if(is.null(cv.fit)) NA else {
#       predicted.age <- predict(cv.fit, as.matrix(df[,sites]), s = 'lambda.min')
#       if(age.transform == 'ln') predicted.age <- exp(predicted.age)
#       median(abs(df$age.best - exp(predict(cv.fit, as.matrix(df[,sites]), s = 'lambda.min'))))
# #      cv.fit$cvm[cv.fit$lambda == cv.fit$lambda.min]
#     } 
#   }, mc.cores = 6) |> 
#     unlist() |> 
#     median(na.rm = TRUE)
# }

glmnet.optim <- do.call(rbind, lapply(c(2, 3, 4), function(minCR){
  do.call(rbind, lapply(weight.schemes, function(weight){
    print(paste0('minCR = ', minCR))
    print(paste0('weight = ', weight))
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
      
      #  find optimal alpha (lowest median MAE at minimum lambda)
      opt.alpha <- optim(
        par = c(alpha = 0.3),
        fn = median.MAE.min,
        method = "Brent",
        lower = 0.001,
        upper = 0.999,
        df = train.df,
        sites = sites,
        age.transform = age.transform,
        nrep = nrep
      )
      return(data_frame(s = sites.2.use, wt = weight, CR = minCR, age.transform = age.transform, a = opt.alpha$par))
    }))
  }))
}))
saveRDS(glmnet.optim, file = 'R/glmnet_tuning/optim.alpha.rds')


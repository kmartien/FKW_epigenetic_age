library(tidyverse)
library(glmnet)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")

optimum.alpha <- readRDS('R/glmnet_tuning/optim.alpha.rds')

nrep <- 1000

age.df$ci.wt <- calc.ci.wt(age.df)

model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>
  left_join(
    rownames_to_column(logit.meth, var = 'swfsc.id')
  ) |> 
  column_to_rownames('swfsc.id') 


# standard 10-fold cv.glmnet run at given alpha
CVglmnet10K <- function(df, sites, age.transform, alpha) {
  if (age.transform == 'ln') df$age.best <- log(df$age.best)
  cv.glmnet(
    x = as.matrix(df[, sites]),
    y = df$age.best,
    weights = df$wt,
    alpha = alpha,
  ) 
}

glmnet.chosen.sites <- 
  lapply(1:nrow(optimum.alpha), function(i){
    print(i)
    print(date())
    minCR <- optimum.alpha$CR[i]
    weight <- optimum.alpha$wt[i]
    age.transform <- optimum.alpha$age.transform[i]
    alpha <- optimum.alpha$a[i]
    model.df <- model.df |> 
      mutate(
        wt = if(weight == 'ci.wt') ci.wt else {
          if (weight == 'CR') age.confidence else {
            if (weight == 'sn.wt') confidence.wt else 1
          }
        })
    train.df <- filter(model.df, age.confidence >= minCR)
    sites <- sites.to.keep
    # multiple LOO runs of the glmnet model at optimum alpha
    fit <- parallel::mclapply(1:nrep, function(i) {
      tryCatch({
        CVglmnet10K(
          train.df, 
          sites.to.keep, 
          age.transform,
          alpha
        )
      }, error = function(e) NULL)
    }, mc.cores = 10)
    site.incl.count <- do.call(bind_rows,
            lapply(fit, function(i){
              coeffs <- coef(i, s = 'lambda.min') 
              data.frame(sites = coeffs@Dimnames[[1]][coeffs@i+1])
            })) |> 
      group_by(sites) |> 
      summarise(count = n()) 
    return(bind_cols(site.incl.count, cr = minCR, wt = weight, age.transform = age.transform, a = alpha))
  }) |> bind_rows()

saveRDS(glmnet.chosen.sites, file = 'R/glmnet_tuning/glmnet.chosen.sites.rds')


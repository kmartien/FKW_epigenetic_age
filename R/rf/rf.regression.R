library(randomForest)
library(rfPermute)
library(dplyr)
library(tidyverse)
source('R/misc_funcs.R')
source('R/tuning_funcs.R')
load("data/age_and_methylation_data.rdata")
load("data/model.params.rda")

nrep <- 1000
ncores <- 10

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

  if(file.exists(file = paste0('R/glmnet_tuning/optim_alpha_', description, '.rds'))){
    opt.alpha <- readRDS(paste0('R/glmnet_tuning/optim_alpha_', description, '.rds'))
  } else {
    param.df  <- rf.param.grid.search(model.df, sites, age, ntree)
    save(param.df, file = paste0('R/rf_tuning/rf_optim_params_', description, '.grid.search.rda'))
  }
  rf.data.age <- select(data.complete, -c(id, wt))
  wt <- data.complete$wt
  
  date()
  rf.age <- rfPermute(
    age.best ~ .,
    rf.data.age,
    ntree = 2000,
    num.rep = 2,
    weights = wt
  )
  date()
  data.complete$predicted.age <- rf.age$rf$predicted
  rf.res <- select(data.complete, c(id, age.best, wt, predicted.age)) %>% 
    left_join(select(all.samples, c(id, sex, age.confidence))) %>% 
    mutate(error = predicted.age - age.best)
  
  mae <- do.call(rbind, lapply(training.min.CR:5, function(cr){
    x <- filter(rf.res, age.confidence == cr)
    age.error <- mutate(x, age.error = abs(error)) %>%
      select(age.error) 
    cor.coeff <- round(cor.test(x$age.best, x$predicted.age, method = "pearson")$estimate,2)
    return(c(CR = cr, MAE = median(age.error$age.error), corr = cor.coeff))
  }))
  
  write.csv(mae, file =paste0("results-raw/rf.mae-", description, ".csv"))
  
  save(rf.age, rf.res, mae, file = paste0("results/rf.regression-", description, ".rda"))
  
  # ImpData <- as.data.frame(importance(rf.age$rf)) %>% filter(`%IncMSE` > 0)
  # ImpData$Var.Names <- row.names(ImpData)
  # 
  # ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  #   geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  #   geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  #   theme_light() +
  #   coord_flip() +
  #   theme(
  #     legend.position="bottom",
  #     panel.grid.major.y = element_blank(),
  #     panel.border = element_blank(),
  #     axis.ticks.y = element_blank()
  #   )
  # 
  # error.hist <- loov.hist(rf.res, min.cr = 2)
  # jpeg(filename = paste0("results-raw/rf.regression.error.histogram", description, ".jpg"))
  # error.hist
  # dev.off()
  # 
  # jpeg(filename = paste0("results-raw/rf.regression.plot.", description, ".jpg"))
  # plot.loov.res(rf.res, min.CR = 2)
  # dev.off()
  # 
  #imp <- rownames_to_column(data.frame(rf.age$rf$importance), var = "site") %>% select(-IncNodePurity)
  #jpeg(filename = paste0("results-raw/rf.regression.site.importance.", description, ".jpg"))
  #plot.site.importance(imp)
  #dev.off()
})


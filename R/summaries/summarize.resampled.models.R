library(tidyverse)
library(ggplot2)
library(gridExtra)
#source("R/plotting.funcs.R")
source('R/misc_funcs.R')
load("data/model.params.rda")
load('data/age_and_methylation_data.rdata')

model.params <- model.params |> 
  rename(sites = site.select.regr.meth) |> 
  mutate(minCR = paste0('minCR', training.cr),
         site.select.cr = paste0('cr', site.select.cr))

# subset all data
age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep) |> 
  mutate(best.age.bin = ifelse(age.best < 10, 0, ifelse(age.best < 25, 10, 25)))


minCRs <- c(2, 3, 4)
sites.2.use <- c('RFsites', 'glmnet', 'Allsites')
wts.2.use <- c('CR', 'ci.wt', 'none')

res.files <- 
  lapply(c('svm', 'gam', 'glmnet', 'rf'), function(m){
    fnames <- c(list.files(paste0('R/', m), pattern = paste0(m, '_ranAge')))
    bind_cols(fnames, do.call(rbind, strsplit(fnames, split = '_')))
  }) |> bind_rows()
names(res.files) <- c('fname', 'method', 'resample', 'minCR', 'sites', 'site.select.cr', 'age.transform', 'weight')
res.files$weight <- substr(res.files$weight, 1, nchar(res.files$weight) - 4)
res.files <- left_join(res.files, select(model.params, -training.cr)) |> 
  relocate(model)


#summarizing predicted age distributions from resampled runs

ran.age.distributions <- 
  left_join(
    res.files,
    lapply(1:nrow(res.files), function(f){
      data.frame(fname = res.files$fname[f], readRDS(paste0('R/', res.files$method[f], '/', res.files$fname[f])))
    }) |> bind_rows(),
    by = 'fname'
  ) |> select(-c(lci, uci, k))

pred.ran <- 
  ran.age.distributions |> 
  filter(weight %in% wts.2.use, sites %in% sites.2.use) |> 
  group_by(model, method, age.transform, sites, weight, minCR, site.select.cr, swfsc.id) |> 
  summarise(age.pred = modeest::venter(age.pred)) |> 
  left_join(select(age.df, c(swfsc.id, age.best))) |> 
  mutate(resid = age.pred - age.best,
         dev = abs(resid),
         sq.err = resid ^ 2,
         best.age.bin = ifelse(age.best < 10, 0, ifelse(age.best < 25, 10, 25)),
         pred.age.bin = ifelse(age.pred < 10, 0, ifelse(age.pred < 25, 10, 25))
  ) 

in.HDI <- ran.age.distributions |>
  group_by(model, method, minCR, resample, sites, site.select.cr, weight, age.transform, swfsc.id) |> 
  summarise(as.data.frame(rbind(HDInterval::hdi(age.pred)))
  ) |> 
  left_join(pred.ran) |> 
  left_join(select(age.df, c(swfsc.id, age.best, age.confidence))) |> 
  filter(age.confidence > 3) |> 
  mutate(in.ci = age.best >= lower & age.best <= upper) |>
  group_by(method, model, resample) |>
  summarise(num.in.ci = sum(in.ci))
in.HDI

MAE.all <- 
  left_join(
    pred.ran |> 
      left_join(age.df) |> 
      filter(age.confidence %in% c(4,5)) |> 
      group_by(model, method, sites, weight, minCR, site.select.cr, age.transform) |> 
      summarise(MAE = round(median(dev), 2),
                lci = round(quantile(dev, probs = c(.25)),2),
                uci = round(quantile(dev, probs = c(.75)),2),
                Corr = round(cor.test(age.best, age.pred, method = "pearson")$estimate,2),
                mean.resid = round(mean(resid), 2)),
    pred.ran |> 
      left_join(age.df) |> 
      filter(age.confidence %in% c(4,5)) |> 
      group_by(model, method, sites, weight, minCR, site.select.cr, age.transform, best.age.bin) |> 
      summarise(MAE = round(median(dev), 2),
                mean.resid = round(mean(resid), 2)) |> 
      pivot_wider(names_from = best.age.bin, values_from = c(MAE, mean.resid))
  )

age.dist.gam.ranAgeMeth <- 
  filter(ran.age.distributions, 
         method == 'svm' & model == 21) |> 
  left_join(select(age.df, c(swfsc.id, age.confidence, age.best)))
composite.age.dist <- lapply(1:nrow(age.df), function(i) {
  crcAgeDist(age.df$swfsc.id[i], age.df, type = 'density') |> 
    mutate(density = ifelse(age < age.df$age.min[i], 0, density),
           density = ifelse(age > age.df$age.max[i], 0, density)) |> 
    filter(age > age.df$age.min[i]-0.1 & age < age.df$age.max[i]+0.1)
}) |>
  bind_rows() |> 
  left_join(select(age.df, c(swfsc.id, age.confidence)))
age.dist <- lapply(2:5, function(cr){ 
  ggplot(filter(age.dist.gam.ranAgeMeth, age.confidence == cr)) +
    geom_histogram(aes(x = age.pred, after_stat(density), fill = factor(age.confidence)), show.legend = FALSE) + 
    geom_vline(aes(xintercept = age.best), linetype = 'dashed', color = 'gray50') +
    geom_line(aes(x = age, y = density), data = filter(composite.age.dist, age.confidence == cr), color = 'gray50') +
#    geom_vline(aes(xintercept = age.pred), data = filter(pred.GAM.ranAgeMeth, age.confidence == cr)) +
    #    geom_segment(aes(x = min(age), xend = min(age), y = 0, yend = age[1]), data = filter(composite.age.dist, age.confidence == cr), linetype = 'dotted', color = 'gray50') +
    #    geom_vline(aes(xintercept = age.max), data = filter(pred.GAM.ranAgeMeth, age.confidence == cr), linetype = 'dotted', color = 'gray50') +
    scale_fill_manual(values = conf.colors) +
    facet_wrap(~swfsc.id, scales = 'free', nrow = 6)
})
pdf(file = 'R/summaries/forPerth/predicted.age.distributions.gam.ranAgeMeth.pdf')#, height = 2000, width = 2000)
age.dist
dev.off()

ggplot(filter(age.dist.gam.ranAgeMeth, swfsc.id == 'z0033906')) +
  geom_histogram(aes(x = age.pred, after_stat(density), fill = factor(age.confidence)), show.legend = FALSE) + 
  geom_vline(aes(xintercept = age.best), linetype = 'dashed', color = 'gray50') +
  geom_line(aes(x = age, y = density), data = filter(composite.age.dist, swfsc.id == 'z0033906'), color = 'gray50') +
  #    geom_vline(aes(xintercept = age.pred), data = filter(pred.GAM.ranAgeMeth, age.confidence == cr)) +
  #    geom_segment(aes(x = min(age), xend = min(age), y = 0, yend = age[1]), data = filter(composite.age.dist, age.confidence == cr), linetype = 'dotted', color = 'gray50') +
  #    geom_vline(aes(xintercept = age.max), data = filter(pred.GAM.ranAgeMeth, age.confidence == cr), linetype = 'dotted', color = 'gray50') +
  scale_fill_manual(values = conf.colors) 
#  facet_wrap(~swfsc.id, scales = 'free')

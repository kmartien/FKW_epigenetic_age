library(tidyverse)
library(ggplot2)
library(gridExtra)
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

methods <- c('svm', 'gam', 'glmnet', 'rf')
minCRs <- c(2, 3, 4)
sites.2.use <- c('RFsites', 'glmnet', 'Allsites')
wts.2.use <- c('CR', 'ci.wt', 'none')

res.files <- 
  lapply(methods, function(m){
  fnames <- c(list.files(paste0('R/', m), pattern = paste0(m, '_best')))
#              list.files(paste0('R/', m), pattern = paste0(m, '_best')))
      bind_cols(fnames, do.call(rbind, strsplit(fnames, split = '_')))
}) |> bind_rows()
names(res.files) <- c('fname', 'method', 'resample', 'minCR', 'sites', 'site.select.cr', 'age.transform', 'weight')
res.files$weight <- substr(res.files$weight, 1, nchar(res.files$weight) - 4)
res.files <- left_join(res.files, select(model.params, -training.cr)) |> 
  relocate(model)

pred <- 
  left_join(
  res.files, 
  lapply(1:nrow(res.files), function(f){
    data.frame(fname = res.files$fname[f], readRDS(paste0('R/', res.files$method[f], '/', res.files$fname[f])))
  }) |> bind_rows(),
  by = 'fname'
) |> 
  filter(weight %in% wts.2.use, sites %in% sites.2.use) |> 
  group_by(model, method, age.transform, sites, weight, minCR, site.select.cr, swfsc.id) |> 
  summarise(age.pred = modeest::venter(age.pred)) |> 
  left_join(select(age.df, c(swfsc.id, age.best))) |> 
  mutate(resid = age.pred - age.best,
         dev = abs(resid),
         sq.err = resid ^ 2,
#         model = paste(minCR, sites, site.select.cr, age.transform, weight, sep = '_'),
         best.age.bin = ifelse(age.best < 10, 0, ifelse(age.best < 25, 10, 25)),
         pred.age.bin = ifelse(age.pred < 10, 0, ifelse(age.pred < 25, 10, 25))
  )
saveRDS(pred, file = 'results/all_models_predictions.rds')

MAE.all <- 
  left_join(
    pred |> 
      left_join(age.df) |> 
      filter(age.confidence %in% c(4,5)) |> 
      group_by(model, method, sites, weight, minCR, site.select.cr, age.transform) |> 
  summarise(MAE = round(median(dev), 2),
                lci = round(quantile(dev, probs = c(.25)),2),
                uci = round(quantile(dev, probs = c(.75)),2),
                Corr = round(cor.test(age.best, age.pred, method = "pearson")$estimate,2),
                mean.resid = round(mean(resid), 2)),
    pred |> 
      left_join(age.df) |> 
      filter(age.confidence %in% c(4,5)) |> 
      group_by(model, method, sites, weight, minCR, site.select.cr, age.transform, best.age.bin) |> 
      summarise(MAE = round(median(dev), 2),
                mean.resid = round(mean(resid), 2)) |> 
      pivot_wider(names_from = best.age.bin, values_from = c(MAE, mean.resid))
  )
write.csv(MAE.all, file = 'results_raw/MAE.all.csv', row.names = FALSE)

top.models <- lapply(methods, function(m){
  temp <- filter(MAE.all, method == m & weight != 'ci.wt') |> 
    arrange(MAE) |> 
    select(model)
  return(temp$model[1:5])
})
names(top.models) <- methods
saveRDS(top.models, file = 'data/top.models.rds')

# Best model plots -----------------------------------------------
# Plot the results of the best model

load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/archive_data.rda")

df <- pred |> 
  filter(model == 21 & method == 'svm') |> 
  left_join(age.df) |> 
  left_join(
    archive_data |> 
      rename(swfsc.id = LABID) |> 
      select(swfsc.id, Collection_Method)
    )

p.loov <- 
  ggplot(filter(df, age.confidence > 3), aes(x = age.best, y = age.pred)) +
  geom_point(aes(col = as.character(age.confidence)), size = 3) +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  labs(x = "Age.best", y = "Predicted age") +
  scale_color_manual(values = conf.colors, name = "Confidence") +
  theme(
    text = element_text(size = 20),
    legend.position = c(0.2, 0.8)
  )
jpeg(file = 'R/summaries/base.model.plot.jpg', width = 600, height = 600)
p.loov
dev.off()

error_by_collection_method <- df |> 
  filter(age.confidence > 3) |> 
  group_by(Collection_Method) |> 
  summarise(MAE = median(dev),
            median.age = median(age.best))


MAE.by.sex <- df |> 
  group_by(sex) |> 
  summarise(MAE = median(dev))


library(tidyverse)
library(ggplot2)
library(ggpubr)
library(viridis)
library(gridExtra)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')

# subset all data
age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep) |> 
  mutate(best.age.bin = ifelse(age.best < 10, 0, ifelse(age.best < 25, 10, 25)))

MAE.all <- read.csv('results_raw/MAE.all.csv')

#############################################################################
# Model design comparisons
#############################################################################

# effect of age transformation
age.plot <- MAE.all |> 
  ungroup() |> 
  filter(weight != 'ci.wt') |> 
  select(c(method, sites, weight, minCR, site.select.cr, age.transform, MAE)) |> 
  mutate(age.transform = recode(age.transform, 
                                'none' = 'None', 
                                'ln' = 'Log')) |>
  group_by(age.transform, method) |> 
  summarise(med.MAE = median(MAE)) |> 
  ggplot(aes(x = age.transform, y = med.MAE, colour = method, group = method)) +
  geom_line(linewidth = 2) +
  geom_point(size = 4) +
  labs(x = 'Age transformation', y = 'Median MAE', colour = 'Training method') +
  theme(text = element_text(size = 40)) +
  scale_color_manual(values = c("svm" = viridis(20)[1], 
                                "gam" = viridis(20)[7], 
                                "glmnet" = viridis(20)[13], 
                                "rf" = viridis(20)[19]),
                     labels = c("svm" = "SVM", 
                                "gam" = "GAM", 
                                "glmnet" = "ENR", 
                                "rf" = "RFR")
                     )
# site selection
site.plot <- MAE.all |> 
  ungroup() |> 
  filter(weight != 'ci.wt') |> 
  select(c(method, sites, weight, minCR, site.select.cr, age.transform, MAE)) |> 
  mutate(sites = recode(sites, 
                                'Allsites' = 'All sites', 
                                'RFsites' = 'RFsites',
                                'glmnet' = 'ENRsites')) |>  
  group_by(sites, method) |> 
  summarise(med.MAE = median(MAE)) |> 
  ggplot(aes(x = sites, y = med.MAE, colour = method, group = method)) +
  geom_line(linewidth = 2) +
  geom_point(size = 4) +
  labs(x = 'CpG site selection', y = 'Median MAE', colour = 'Training method') +
  theme(text = element_text(size = 40)) +
  scale_color_manual(values = c("svm" = viridis(20)[1], 
                                "gam" = viridis(20)[7], 
                                "glmnet" = viridis(20)[13], 
                                "rf" = viridis(20)[19]),
                     labels = c("svm" = "SVM", 
                                "gam" = "GAM", 
                                "glmnet" = "ENR", 
                                "rf" = "RFR")
  )
# minCR
minCR.plot <- MAE.all |> 
  ungroup() |> 
  filter(weight != 'ci.wt') |> 
  select(c(method, sites, weight, minCR, site.select.cr, age.transform, MAE)) |> 
  mutate(minCR = recode(minCR, 
                                'minCR2' = 'All', 
                                'minCR3' = 'CR3+',
                                'minCR4' = 'CR4+')) |>  
  group_by(minCR, method) |> 
  summarise(med.MAE = median(MAE)) |> 
  ggplot(aes(x = minCR, y = med.MAE, colour = method, group = method)) +
  geom_line(linewidth = 2) +
  geom_point(size = 4) +
  labs(x = 'Training sample set', y = 'Median MAE', colour = 'Training method') +
  theme(text = element_text(size = 40)) +
  scale_color_manual(values = c("svm" = viridis(20)[1], 
                                "gam" = viridis(20)[7], 
                                "glmnet" = viridis(20)[13], 
                                "rf" = viridis(20)[19]),
                     labels = c("svm" = "SVM", 
                                "gam" = "GAM", 
                                "glmnet" = "ENR", 
                                "rf" = "RFR")
  )
# weight
weight.plot <- MAE.all |> 
  ungroup() |> 
  filter(weight != 'ci.wt') |> 
  select(c(method, sites, weight, minCR, site.select.cr, age.transform, MAE)) |> 
  mutate(weight = recode(weight, 
                                'CR' = 'CR', 
                                'none' = 'None')) |>  
  group_by(weight, method) |> 
  summarise(med.MAE = median(MAE)) |> 
  ggplot(aes(x = weight, y = med.MAE, colour = method, group = method)) +
  geom_line(linewidth = 2) +
  geom_point(size = 4) +
  labs(x = 'Sample weighting', y = 'Median MAE', colour = 'Training method') +
  theme(text = element_text(size = 40)) +
  scale_color_manual(values = c("svm" = viridis(20)[1], 
                                "gam" = viridis(20)[7], 
                                "glmnet" = viridis(20)[13], 
                                "rf" = viridis(20)[19]),
                     labels = c("svm" = "SVM", 
                                "gam" = "GAM", 
                                "glmnet" = "ENR", 
                                "rf" = "RFR")
  )

jpeg('results_raw/model_design_comparisons.jpg', width = 1500, height = 1500)
ggarrange(age.plot, site.plot, minCR.plot, weight.plot, 
          nrow = 2,
          ncol = 2,
          labels = c('A', 'B', 'C', 'D'), 
#          font_label = list(size = 40, color = 'red'), 
          common.legend = TRUE, 
          legend = 'bottom')
dev.off()

# plot histogram of residuals in each age group for all models trained by a given method
MAE.all |> 
  select(c(model, method, mean.resid_0, mean.resid_10, mean.resid_25)) |> 
  pivot_longer(cols = c(mean.resid_0, mean.resid_10, mean.resid_25), names_to = 'age.class', values_to = 'mean.resid') |> 
  mutate(age.class = recode(age.class,
                            'mean.resid_0' = '0 - 9',
                            'mean.resid_10' = '10 - 24',
                            'mean.resid_25' = '25+')) |> 
  ggplot() +
  geom_histogram(aes(x = mean.resid)) +
  facet_wrap(method~age.class, nrow = 4)

# plot histogram of MAE in each age group for all models trained by a given method
MAE.all |> 
  select(c(model, method, MAE_0, MAE_10, MAE_25)) |> 
  pivot_longer(cols = c(MAE_0, MAE_10, MAE_25), names_to = 'age.class', values_to = 'MAE') |> 
  mutate(age.class = recode(age.class,
                            'MAE_0' = '0 - 9',
                            'MAE_10' = '10 - 24',
                            'MAE_25' = '25+')) |> 
  ggplot() +
  geom_histogram(aes(x = MAE)) +
  facet_wrap(method~age.class, nrow = 4)

#############################################################################
# Best model results
#############################################################################

load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/archive_data.rda")
pred_ran <- readRDS('results/best_model_pred_ran.rds')
pred <- readRDS('results/all_models_predictions.rds')

df <- pred |> 
  filter(model == 21 & method == 'svm') |> 
  left_join(pred_ran) |> 
  ungroup() |> 
  select(-c(age.transform, sites, weight, minCR, site.select.cr)) |> 
  left_join(age.df) |> 
  left_join(
    archive_data |> 
      rename(swfsc.id = LABID) |> 
      select(swfsc.id, Collection_Method)
  )

# Scatterplot of the best model
p.loov <- df |> 
  mutate(age.best = jitter(age.best, amount = 0.5),
         age.pred = jitter(age.pred, amount = 0.5)) |> 
  filter(age.confidence >3 ) |> 
  ggplot(aes(x = age.best, y = age.pred)) +
  geom_point(aes(col = as.character(age.confidence)), size = 3) +
  geom_errorbarh(aes(xmin = age.min, xmax = age.max, col = as.character(age.confidence)), height = 0) +
  geom_errorbar(aes(ymin = lower, ymax = upper, col = as.character(age.confidence)), width = 0) +
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

# looking for differences in residuals between stranded, biopsied, and tagged samples
error_by_collection_method <- df |> 
  filter(age.confidence > 3) |> 
  group_by(Collection_Method) |> 
  summarise(MAE = median(dev),
            median.age = median(age.best))

# check for differences in MAE by sex
MAE.by.sex <- df |> 
  group_by(sex) |> 
  summarise(MAE = median(dev),
            mean.res = mean(resid))


ran_age_dist <- readRDS('results/best_model_ran_age_distribution.rds')

age.dist.gam.ranAgeMeth <- 
  ran_age_dist |> 
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


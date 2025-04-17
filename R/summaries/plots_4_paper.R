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
  annotation_custom(text_grob(label = paste0(
    'Correlation\nCoefficient = ',
    round(
      cor.test(
      formula = ~ age.best + age.pred, 
      data = filter(df, age.confidence > 3),
      method = "pearson")$estimate,
      2)
  ), size = 20), 
  xmin = 20, ymax = 25,) +
#  cor.coeff <- round(cor.test(dat$age.best, dat$age.pred, method = "pearson")$estimate,2)
  labs(x = "Age.best", y = "Predicted age") +
  scale_color_manual(values = conf.colors, name = "Confidence") +
  theme(
    text = element_text(size = 20),
    legend.position = c(0.2, 0.8)
  )
jpeg(file = 'results_raw/base.model.plot.jpg', width = 600, height = 600)
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

# Compare duplicates

df$date.biopsy <- as.POSIXct(age.df$date.biopsy, format = '%Y-%m-%d %H:%M:%S') |> 
  as.Date()

pair.sum <- lapply(
  select(df, c(crc.id, swfsc.id, date.biopsy)) |> filter(n() > 1, .by = crc.id) |> 
    pull(crc.id) |> 
    unique(),
  function(i){
    inds <- filter(df, crc.id == i) |> arrange(date.biopsy)
    data.frame(
      crc.id = i,
      old.ind = inds$swfsc.id[2], 
      young.ind = inds$swfsc.id[1], 
      old.cr = inds$age.confidence[2],
      young.cr = inds$age.confidence[1],
      actual.age.diff = difftime(inds$date.biopsy[2], inds$date.biopsy[1], units = "days") / 365,
      predicted.age.diff = (inds$age.pred[2] - inds$age.pred[1])
    )
  }) |> bind_rows()

ran_age_dist$iter <- do.call(c, lapply(1:1000, function(i){rep(i, 89)}))

# Calculate difference in predicted ages in each iteration
ran_age_diff <- pair.sum %>%
  mutate(across(c(old.ind, young.ind), as.character)) %>%
  rowwise() %>%
  reframe(
    crc.id = crc.id,
    old.ind = old.ind,
    young.ind = young.ind,
    # Create rows for each iteration
    iter = unique(ran_age_dist$iter),
    age_pred_diff = 
      ran_age_dist %>%
      filter(swfsc.id == old.ind, iter == .data$iter) %>%
      pull(age.pred) -
      ran_age_dist %>%
      filter(swfsc.id == young.ind, iter == .data$iter) %>%
      pull(age.pred)
  ) %>%
  # Select and order desired columns
  select(crc.id, old.ind, young.ind, age_pred_diff, iter) |> 
  left_join(select(pair.sum, crc.id, actual.age.diff))

prop_gt_zero <-
  ran_age_diff |> 
  filter(age_pred_diff > 0) |> 
  group_by(crc.id) |> 
  summarise(positive_diff = n()/1000)

# plot distribution of age differences across iterations for each pair
pair.diff.plot <- 
  ggplot(ran_age_diff) +
  geom_histogram(aes(x = age_pred_diff), position = 'identity') +
  geom_vline(aes(xintercept = actual.age.diff), linewidth = 2) +
  geom_vline(aes(xintercept = 0), linewidth = 2, linetype = 2) +
  geom_text(
    data = prop_gt_zero, 
    aes(
      x = -17, 
      y = 200, 
      label = round(positive_diff, 2)),
    hjust = 0,
    vjust = 0,
    size = 10
  ) +
  labs(x = 'Predicted age difference', y = 'Count') +
  theme(text = element_text(size = 35)) +
  facet_wrap(~crc.id, ncol = 3) +
  theme(
    text = element_text(size = 35),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

jpeg(filename = 'results_raw/pair_diff_plot.jpg', width = 1000, height = 1500)
pair.diff.plot
dev.off()

# compare the predicted age probability distributions of young and old individuals from each pair
pair.plot <- lapply(1:nrow(pair.sum), function(i){
  ids <- pair.sum[i, c('old.ind', 'young.ind')]
  ggplot(filter(age.dist.gam.ranAgeMeth, swfsc.id == ids$old.ind | swfsc.id == ids$young.ind)) +
  geom_histogram(
    aes(x = age.pred, 
        fill = factor(age.confidence),
        alpha = factor(swfsc.id),
        color = factor(swfsc.id)),
#        color = factor(swfsc.id)),
    position = 'identity',
    #      alpha = 0.8,
    show.legend = FALSE) + 
  scale_fill_manual(values = conf.colors) +
#  scale_fill_manual(values = c('gray20', 'white')) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_color_manual(values = c('transparent', 'black')) +
#  scale_color_manual(values = c('transparent', 'darkgray')) +
  geom_vline(
    aes(xintercept = age.pred), 
    data = filter(df, swfsc.id == ids$old.ind),
    color = 'gray60') +
  geom_vline(
    aes(xintercept = age.pred), 
    data = filter(df, swfsc.id == ids$young.ind),
    color = 'gray20') +
  annotation_custom(tableGrob(data.frame(
    Difference = c(as.character(round(pair.sum$actual.age.diff[i], 2)), as.character(round(pair.sum$predicted.age.diff[i],2))),
    row.names = c('Actual', 'Predicted')),
    theme = ttheme_minimal(base_size = 16)), xmax = Inf, ymax = Inf) +
    labs(x = element_blank(), y = element_blank()) +
#  labs(x = 'Predicted age', y = 'Density') +
#  theme_minimal() +
  theme(text = element_text(size = 25))
})
pair.plot$ncol <- 3
pair.plot$bottom <- 'Predicted age'
pair.plot$left <- 'Density'
jpeg(filename = 'results_raw/pair.plot.jpg', width = 1000, height = 1500)
do.call(grid.arrange, pair.plot)
dev.off()  


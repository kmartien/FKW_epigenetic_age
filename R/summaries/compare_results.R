library(tidyverse)
library(ggplot2)
library(gridExtra)
source("R/plotting.funcs.R")
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

write.csv(MAE.all, file = 'results_raw/MAE.all.csv')

# Base model box plot -----------------------------------------------
# Plot the results of the base model for each method

df <- pred |> 
  filter((model == 21 & method == 'svm') | (model == 75 & method == 'gam')) |> 
  left_join(age.df) |> 
  mutate(#mode = factor(model, levels = best.models, ordered = TRUE),
         method = factor(method, levels = c('gam', 'svm')))
plot.titles <- data.frame(label = c('gam', 'svm'), title = c('GAM', 'SVM'))
base.plot <- lapply(c('gam', 'svm'), function(m){
  p <- plot.loov.res(filter(df, method == m), min.CR = 4)$p.loov + 
    ggtitle(plot.titles$title[which(plot.titles$label == m)]) 
  return(p)
})
names(base.plot) <- c('GAM', 'SVM')
base.plot$nrow = 2
jpeg(file = 'R/summaries/base.model.plot.jpg', width = 600, height = 1200)
do.call(grid.arrange, base.plot)
dev.off()

residuals.base.plot <- lapply(c('gam', 'svm', 'glmnet', 'rf'), function(m){
  p <- plot.residuals(filter(df, method == m), min.CR = 4) + 
    ggtitle(plot.titles$title[which(plot.titles$label == m)]) 
  return(p)
})
residuals.base.plot$nrow = 2
jpeg(file = 'R/summaries/residuals.plot.jpg', width = 1200, height = 800)
do.call(grid.arrange, residuals.base.plot)
dev.off()

# residuals.plot <-
#   filter(df, age.confidence > 3) |> 
#   ggplot(aes(x = age.best, y = resid)) +
#   geom_point(aes(color = factor(age.confidence)), size = 3, show.legend = FALSE) +
#   stat_smooth(aes(color = 'black', fill = NULL), 
#               show.legend = FALSE,
#               method = "lm", 
#               formula = y ~ x, 
#               geom = "smooth") +
#   scale_color_manual(values = conf.colors) +
#   geom_abline(slope = 0, color = "black", linewidth = 0.5, linetype = 2) +
#   annotation_custom(tableGrob(data.frame(
#     filter(df, age.confidence > 3) |> summarise(MAE = round(median(dev), 2)) |> pull(MAE),
#     Corr = filter(df, age.confidence > 3) |> summarise(round(cor.test(age.best, age.pred, method = "pearson")$estimate,2))),
#     theme = ttheme_minimal(base_size = 16)), xmax = Inf, ymax = Inf) +
#   labs(x = 'Agebest',
#        y = 'Residual') +  
#   theme_minimal() +
#   theme(text = element_text(size = 20)) +
#   facet_wrap(~method, nrow = 4,
#              labeller = labeller(method = c(
#                'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
#              )))
# jpeg(file = 'R/summaries/residuals.plot.jpg', width = 600, height = 1200)
# residuals.base.plot
# dev.off()
df |> filter(age.best > 24) |> group_by(method) |> filter(age.confidence > 3) |> summarise(mean.resid = mean(resid))

# box.plot <-
#   pred |>
#   filter(model == 'best_minCR4_RFsites_ln_none') |>
#   left_join(age.df) |>
#   filter(age.confidence %in% c(4,5)) |>
#   ggplot() +
#   geom_boxplot(aes(x = method, y = resid)) +
#   ylim(-25, 25) +
#   scale_x_discrete(labels = c('GAM', 'ENR', 'RF', 'SVM')) +
#   labs(x = "Method",
#        y = "Absolute age error (yrs)") +
#   theme(text = element_text(size = 24))
# jpeg(file = 'R/summaries/base.model.boxplot.jpg', width = 700, height = 1200)
# box.plot
# dev.off()

# Site selection box plot -----------------------------------------------

box.plot <- 
  pred |> 
  left_join(age.df) |> 
  filter(model %in% c('best_minCR4_RFsites_ln_none', 
                      'best_minCR4_Allsites_ln_none',
                      'best_minCR4_glmnet.5_ln_none'),
         age.confidence %in% c(4,5),
         method == 'gam') |> 
  mutate(sites = factor(sites, levels = c('RFsites', 'glmnet.5', 'Allsites'), ordered = TRUE)) |> 
  ggplot() +
  geom_boxplot(aes(x = sites, y = resid)) +
  scale_x_discrete(labels = c('RF sites', 'ENR sites', 'All sites')) +
  #  ylim(c(0,20)) +
  labs(x = "CpG set used",
       y = "Absolute age error (yrs)") +  
#  theme_minimal() +
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust=1))
  # facet_wrap(~method, nrow = 1, labeller = labeller(method = c(
  #   'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  # )))
jpeg(file = 'R/summaries/site.selection.boxplot.jpg', width = 400, height = 700)
box.plot
dev.off()


# minCR and weight box plot -----------------------------------------------

box.plot <- 
  pred |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5),
         method == 'gam',
         resample == 'best') |> 
         mutate(weight = factor(weight, levels = c('none', 'CR', 'ci.wt'), ordered = TRUE)) |> 
  ggplot() +
  geom_boxplot(aes(x = weight, y = resid)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_x_discrete(labels = c('Unweighted', 'CR', 'RW')) +
  labs(x = "Weight scheme",
       y = "Residual (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~minCR, nrow = 1, labeller = labeller(method = c(
    'minCR2' = 'All samples', 'minCR3' = 'CR3+ samples', 'minCR4' = 'CR4+ samples')))
jpeg(file = 'R/summaries/minCR.weight.boxplot.jpg', width = 700, height = 700)
box.plot
dev.off()

# summarise RanAge and RanAgeMeth models -----------------------------------------------

box.plot <- 
  pred |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5),
         method %in% c('gam', 'svm'),
         model %in% c('best_minCR4_RFsites_ln_none',
                      'ranAge_minCR4_RFsites_ln_none',
                      'ranAge_minCR3_RFsites_ln_none',
                      'ranAge_minCR2_RFsites_ln_none',
                      'ranAgeMeth_minCR4_RFsites_ln_none',
                      'ranAgeMeth_minCR3_RFsites_ln_none',
                      'ranAgeMeth_minCR2_RFsites_ln_none')) |>
  mutate(model = factor(model, levels = 
                          c('best_minCR4_RFsites_ln_none',
                            'ranAge_minCR4_RFsites_ln_none',
                            'ranAge_minCR3_RFsites_ln_none',
                            'ranAge_minCR2_RFsites_ln_none',
                            'ranAgeMeth_minCR4_RFsites_ln_none',
                            'ranAgeMeth_minCR3_RFsites_ln_none',
                            'ranAgeMeth_minCR2_RFsites_ln_none'),
                        ordered = TRUE)) |> 
  ggplot() +
  geom_boxplot(aes(x = model, y = resid)) +
  #scale_x_discrete(labels = c('None', 'HCsamps, Age', 'Allsamps, Age', 'HCsamps, Age&Meth', 'Allsamps, Age&Meth')) +
#  ylim(c(0,25)) +
  labs(x = "Resampling",
       y = "Absolute age error (yrs)") + 
  theme_minimal() +
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~method, nrow = 1, labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/resampling.boxplot.jpg', width = 700, height = 900)
box.plot
dev.off()

# Scatterplot and residuals plot for GAM and SVM ranAgeMeth
df <- 
  pred |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5),
         method %in% c('gam', 'svm'),
         model %in% c('best_minCR4_RFsites_ln_none',
                      'ranAgeMeth_minCR4_RFsites_ln_none')) |> 
  mutate(age.confidence = factor(age.confidence))

ranAgeMeth.plot <- ggplot(df, aes(x = age.best, y = age.pred)) +
  geom_point(aes(col = age.confidence), size = 3, show.legend = FALSE) +
  stat_smooth(aes(color = 'black', fill = NULL), 
              show.legend = FALSE,
              method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  scale_color_manual(values = conf.palette, name = "Confidence") +
  labs(x = "Age.best", y = "Predicted age") +
  xlim(0,40) + ylim(0,58) +
  theme_minimal() +
  theme(
    text = element_text(size = 20)
  ) +
  facet_wrap(resample~method, nrow = 4)
jpeg(file = 'R/summaries/ranAgeMeth.minCR.scatterplot.jpg', width = 600, height = 1200)
ranAgeMeth.plot
dev.off()

residuals.plot <-
  ranAgeMeth.plot <- ggplot(df, aes(x = age.best, y = resid)) +
  geom_point(aes(col = age.confidence), size = 3, show.legend = FALSE) +
  stat_smooth(aes(color = 'black', fill = NULL), 
              show.legend = FALSE,
              method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  geom_abline(slope = 0, color = "black", linewidth = 0.5, linetype = 2) +
  scale_color_manual(values = conf.palette, name = "Confidence") +
  labs(x = "Age.best", y = "Residual minCR4") +
  ylim(-20,20) +
  theme_minimal() +
  theme(
    text = element_text(size = 20)
  ) +
  facet_wrap(resample~method, nrow = 4)
jpeg(file = 'R/summaries/ranAgeMeth.residual.jpg', width = 600, height = 1200)
residuals.plot
dev.off()


#summarizing predicted age distributions from resampled runs

ran.age.distributions <- 
  left_join(
    filter(res.files, resample != 'best',
           method %in% c('gam', 'svm')),
    lapply(1:nrow(res.files), function(f){
      data.frame(fname = res.files$fname[f], readRDS(paste0('R/', res.files$method[f], '/', res.files$fname[f])))
    }) |> bind_rows(),
    by = 'fname'
  ) |> select(-c(lci, uci, k)) |>
  filter(weight == 'none', sites == 'RFsites') 

pred.ran <- ran.age.distributions |>
  group_by(method, minCR, resample, sites, weight, swfsc.id) |> 
  summarise(as.data.frame(rbind(HDInterval::hdi(age.pred)))
  ) |> 
  left_join(pred) 

in.HDI <- pred.ran |> 
  left_join(select(age.df, c(swfsc.id, age.best, age.confidence))) |> 
  filter(age.confidence > 3) |> 
  mutate(in.ci = age.best >= lower & age.best <= upper) |>
  group_by(method, minCR, resample) |>
  summarise(num.in.ci = sum(in.ci))
in.HDI

# Summarizing GAM ranAgeMeth model

# MAE by PREDICTED age for ranAgeMeth
breaks <- c(0,10,25,40)
pred.GAM.ranAgeMeth <- 
  pred.ran |> 
  left_join(age.df) |> 
  filter(method == 'svm',
         model == 21) |> 
  mutate(best.age.bin = ifelse(age.best < 10, 0, ifelse(age.best < 25, 10, 25)),
         pred.age.bin = ifelse(age.pred < 10, 0, ifelse(age.pred < 25, 10, 25)))
MAE.by.age <- 
  left_join(
    MAE.pred <- 
      filter(pred.GAM.ranAgeMeth, age.confidence > 3) |> 
      group_by(pred.age.bin) |> 
      summarise(MAE.pred = round(median(dev), 2)) |> 
      rename(age.bin = pred.age.bin),
    MAE.best <- filter(pred.GAM.ranAgeMeth, age.confidence > 3) |> 
      group_by(best.age.bin) |> 
      summarise(MAE.best = round(median(dev), 2)) |> 
      rename(age.bin = best.age.bin),
    by = 'age.bin'
  ) 
write.csv(MAE.by.age, file = paste0("R/summaries/MAE.by.age.gam.ranAgeMeth.csv"))

# scatterplot with error bars
GAM.ranAgeMeth.plot <- pred.GAM.ranAgeMeth |> 
#  filter(age.confidence > 3) |> 
  mutate(age.best = jitter(age.best),
         age.pred = jitter(age.pred),
    age.confidence = factor(age.confidence)) |> 
  ggplot(aes(x = age.best, y = age.pred)) +
  geom_point(aes(col = age.confidence), size = 3, show.legend = FALSE, position = 'jitter') +
  geom_segment(aes(x = age.best, xend = age.best, y = lower, yend = upper, color = age.confidence), alpha = 0.5, position = 'jitter') +
  geom_segment(aes(x = age.min, xend = age.max, y = age.pred, yend = age.pred, color = age.confidence), alpha = 0.5, position = 'jitter') +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  scale_color_manual(values = conf.colors, name = "Confidence") +
  labs(x = "Age.best", y = "Predicted age") +
#  xlim(0,40) + ylim(0,58) +
  facet_wrap(~age.confidence, nrow = 2) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    legend.position = 'top',
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
jpeg(file = 'R/summaries/GAM.ranAgeMeth.scatterplot.jpg', width = 600, height = 1200)
GAM.ranAgeMeth.plot
dev.off()

age.dist.gam.ranAgeMeth <- 
  filter(ran.age.distributions, 
  method == 'gam' & minCR == 'minCR4' & resample == 'ranAgeMeth') |> 
  left_join(select(age.df, c(swfsc.id, age.confidence)))
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
    geom_line(aes(x = age, y = density), data = filter(composite.age.dist, age.confidence == cr), color = 'gray50') +
    geom_vline(aes(xintercept = age.pred), data = filter(pred.GAM.ranAgeMeth, age.confidence == cr)) +
    geom_vline(aes(xintercept = age.best), data = filter(pred.GAM.ranAgeMeth, age.confidence == cr), linetype = 'dashed', color = 'gray50') +
#    geom_segment(aes(x = min(age), xend = min(age), y = 0, yend = age[1]), data = filter(composite.age.dist, age.confidence == cr), linetype = 'dotted', color = 'gray50') +
#    geom_vline(aes(xintercept = age.max), data = filter(pred.GAM.ranAgeMeth, age.confidence == cr), linetype = 'dotted', color = 'gray50') +
    scale_fill_manual(values = conf.colors) +
  facet_wrap(~swfsc.id, scales = 'free')
})
pdf(file = 'R/summaries/predicted.age.distributions.gam.ranAgeMeth.pdf')#, height = 2000, width = 2000)
age.dist
dev.off()


# compare duplicates -----------------------------------------------

age.df$date.biopsy <- as.POSIXct(age.df$date.biopsy, format = '%Y-%m-%d %H:%M:%S') %>% 
  as.Date()

pair.sum <- lapply(
  select(age.df, c(crc.id, swfsc.id, date.biopsy)) |> 
    filter(n() > 1, .by = crc.id) |> 
    pull(crc.id) |> 
    unique(), 
  function(i){
    inds <- filter(pred.GAM.ranAgeMeth, crc.id == i) %>% arrange(date.biopsy)
    data.frame(
      crc.id = i,
      old.ind = inds$swfsc.id[2], 
      young.ind = inds$swfsc.id[1], 
      old.cr = inds$age.confidence[2],
      young.cr = inds$age.confidence[1],
      actual.age.diff = difftime(inds$date.biopsy[2], inds$date.biopsy[1], units = "days")/365,
      predicted.age.diff = inds$age.pred[2] - inds$age.pred[1]
    )
  }) |> bind_rows() 

pair.plot <- lapply(1:nrow(pair.sum), function(i){
  ids <- pair.sum[i, c('old.ind', 'young.ind')]
  ggplot(filter(age.dist.gam.ranAgeMeth, swfsc.id == ids$old.ind | swfsc.id == ids$young.ind)) +
    geom_histogram(
      aes(x = age.pred, 
          fill = factor(age.confidence),
          alpha = factor(swfsc.id),
          color = factor(swfsc.id)),
      position = 'identity',
#      alpha = 0.8,
      show.legend = FALSE) + 
#    scale_fill_manual(values = c('gray20', 'gray60')) +
    scale_fill_manual(values = conf.colors) +
    scale_alpha_manual(values = c(1, 0.8)) +
    scale_color_manual(values = c('transparent', 'gray')) +
    geom_vline(
      aes(xintercept = age.pred), 
      data = filter(pred.GAM.ranAgeMeth, swfsc.id == ids$old.ind),
      color = 'gray60') +
    geom_vline(
      aes(xintercept = age.pred), 
      data = filter(pred.GAM.ranAgeMeth, swfsc.id == ids$young.ind),
      color = 'gray20') +
    annotation_custom(tableGrob(data.frame(
      Difference = c(as.character(round(pair.sum$actual.age.diff[i], 2)), as.character(round(pair.sum$predicted.age.diff[i],2))), 
      row.names = c('Actual', 'Predicted')),
      theme = ttheme_minimal(base_size = 16)), xmax = Inf, ymax = Inf) +
    labs(x = 'Predicted age', y = 'Density') +
    theme_minimal() +
    theme(text = element_text(size = 16))
})
pair.plot$ncol <- 4
jpeg(filename = 'R/summaries/pair.plot.jpg', width = 2000, height = 1500)
do.call(grid.arrange, pair.plot)
dev.off()

# MAE by age and CR -----------------------------------------------
breaks <- c(0,10,25,40)
df <- pred |>
  left_join(age.df) |> 
  filter(model == 'best_minCR4_RFsites_ln_none' & age.confidence > 3)

MAE.by.age <-  do.call(rbind, lapply(1:(length(breaks)-1), function(a){
  filter(df, age.best >= breaks[a]) %>% 
    filter(age.best
           < breaks[a+1]) %>% 
    group_by(method) |> 
    summarise(MAE = median(dev)) |> 
    mutate(age_bin = breaks[a]) #%>% 
})) |> 
  pivot_wider(names_from = method, values_from = MAE)
write.csv(MAE.by.age, file = paste0("R/summaries/MAE.by.age.base.models.csv"))

# Check base GAM results ------------------------------------------------

base.gam <- readRDS(file = 'R/gam/gam_best_minCR4_RFsites_ln_none.rds') |> 
  left_join(age.df) |> 
  mutate(in.ci = ifelse(age.best >= lci & age.best <= uci, 1, 0))

filter(base.gam, age.confidence > 3) |> 
  group_by(age.confidence) |> 
  summarise(sum(in.ci))

###########################################################################
# code I'm not currently using 


box.plot <- 
  pred.2.plot |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  mutate(age.conficence = as.factor(age.confidence)) |> 
  mutate(minCR = as.factor(minCR)) |> 
  ggplot() +
  geom_boxplot(aes(x = model, y = dev)) +
  labs(x = "",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_blank()) +
  facet_wrap(~method, nrow = 2,labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/model.boxplot.jpg', width = 1200, height = 1800)
box.plot
dev.off()


# plot age.best vs. inv.var
plots <-
  lapply(c('CR', 'sn.wt', 'inv.var'), function(weight){
    age.df |> 
      mutate(wt = if (weight == 'CR') age.confidence else {if(weight == 'inv.var') 1/age.var else confidence.wt}) |> 
#      mutate(age.confidence = factor(age.confidence)) |> 
      ggplot() +
      geom_point(aes(x = age.best, y = wt))#, colour = age.confidence)) +
#      scale_fill_manual(values = conf.colors)
  })



base.res <-   pred |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  filter(weight == 'none' & minCR == 4)
  
MAE.base <- 
  base.res |> 
  group_by(method, sites, weight, minCR) |> 
  summarise(MAE = median(dev))

box.plot <- 
  base.res |> 
  mutate(age.conficence = as.factor(age.confidence)) |> 
  mutate(minCR = as.factor(minCR)) |> 
  ggplot() +
  geom_boxplot(aes(x = sites, y = dev)) +
  labs(x = "CpG site selection",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, hjust=1)) +
#  facet_wrap(~method, nrow = 1)
  facet_wrap(~method, nrow = 1,labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/site.selection.jpg', width = 900, height = 1200)
box.plot
dev.off()

# sample selection and weighting (sites = RFsites)

minCR.res <-   pred |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  filter(sites == 'RFsites') |> 
  filter(minCR == 2 | weight == 'none' & minCR == 4)

MAE.minCR <- 
  minCR.res |> 
  group_by(method, sites, weight, minCR) |> 
  summarise(MAE = median(dev))

box.plot <- 
  minCR.res |> 
  mutate(age.conficence = as.factor(age.confidence)) |> 
  mutate(minCR = as.factor(minCR)) |> 
  mutate(minCR.wt = paste0('minCR', minCR, '_', weight)) |> 
  ggplot() +
  geom_boxplot(aes(x = minCR.wt, y = dev)) +
  labs(x = "minCR",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, hjust=1)) +
#    facet_wrap(~method, nrow = 1)
  facet_wrap(~method, nrow = 1, labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/minCR.jpg', width = 1800, height = 900)
box.plot
dev.off()

# histogram of MAE by method ------------------------------------------------
MAE.hist <- MAE.all |> 
  ggplot() +
  geom_histogram(aes(x = MAE, fill = method)) 
MAE.hist
-----------------------------------------------------------------------------

plots <- lapply(1:length(dat), function(i){
  p <- plot.loov.res(dat[[i]], min.CR = 4)
  p$p.loov$labels$title <- paste0(names(dat)[i])
  return(p$p.loov)
})

#plots[[1]]$labels$title <- "ENR optimized alpha = 0.1"
#plots[[2]]$labels$title <- "ENR alpha = 0.5"
plots$nrow <- 2
jpeg(file = paste0("R/summaries/plainjane.regression.plots--minCR", minCR, "_", sites.2.use, ".jpg"), width = 1000, height = 800)
do.call(grid.arrange, plots)
dev.off()

# box plots
age.errors.long <- do.call(bind_rows, lapply(1:length(dat), function(i){
  bind_cols(method = names(dat)[i], dat[[i]])
}))
age.errors.long$CR <- as.factor(age.errors.long$age.confidence)
box.plot <- 
  ggplot(age.errors.long) +
  geom_boxplot(aes(x = CR, y = dev)) +
  labs(x = "Confidence Rating",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, hjust=1)) +
  #  facet_wrap(~training.cr, nrow = 3, labeller = labeller(training.cr = training.cr.labs))  
  facet_wrap(~method, nrow = 2)  

# MAE for CR=5 males vs. females
MAE.by.sex <- do.call(rbind, lapply(1:length(dat), function(i){
  filter(dat[[i]], age.confidence == 5) %>% group_by(sex) %>% 
    summarise(MAE = median(dev)) %>% bind_cols(method = names(dat)[i])
}))

# MAE by age and CR
breaks <- c(0,10,25,40)
MAE.by.age <- do.call(cbind, lapply(1:length(dat), function(i){
  do.call(rbind, lapply(1:(length(breaks)-1), function(a){
    filter(dat[[i]], age.best >= breaks[a]) %>% filter(age.best < breaks[a+1]) %>% 
      # group_by(age.confidence) %>% #use this line to compare across CRs
      filter(age.confidence >= 4) %>% #use this line to combined CR4and5
      summarise(MAE = median(dev)) %>% bind_cols(age_bins = breaks[a])
  }))
}))
write.csv(MAE.by.age, file = paste0("R/summaries/MAE.by.age.across.methods-minCR", minCR, "_", sites.2.use, ".csv"))

#compare duplicates  
dupe.sum <- do.call(bind_rows, lapply(1:length(dat), function(m){
  dupes <- left_join(dat[[m]], select(age.df, c(crc.id, swfsc.id, date.biopsy))) %>% 
    filter(n() > 1, .by = crc.id)
  do.call(rbind, lapply(unique(dupes$crc.id), function(i){
    inds <- filter(dupes, crc.id == i) %>% arrange(age.best)
    actual.diff <- difftime(inds$date.biopsy[2], inds$date.biopsy[1], units = "days")/365 %>% 
      as.numeric()
    predicted.diff <- inds$age.pred[2] - inds$age.pred[1]
    return(data.frame(method = names(dat)[m], crc.id = i, actual.diff = actual.diff, predicted.diff = predicted.diff, age.confidence = as.character(inds$age.confidence[1])))
  }))
}))
pair.plot <- ggplot(dupe.sum) +
  geom_point(aes(x = actual.diff, y = predicted.diff, colour = age.confidence), size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_colour_manual(values = conf.colors) +
  labs(x = "Actual age difference", y = "Predicted age difference") +
  theme(text = element_text(size = 24)) +
  facet_wrap(~method, nrow = 2)
jpeg(filename = paste0("R/summaries/pair.plot_minCR", minCR, "_", sites.2.use, ".jpg"), width = 960, height = 960)
pair.plot
dev.off()


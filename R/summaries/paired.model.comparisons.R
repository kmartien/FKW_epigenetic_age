library(tidyverse)
library(ggplot2)
library(ggpubr)
library(viridis)
library(gridExtra)

MAE.all <- read.csv('results_raw/MAE.all.csv')

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

# Rf to determine which model choice parameters are important
temp <- MAE.all |> 
  ungroup() |> 
  select(c(method, sites, weight, minCR, site.select.cr, age.transform, MAE)) |> 
  mutate(method.num = ifelse(method == 'gam', 1, 
                             ifelse(method == 'svm', 2, 
                                    ifelse(method == 'glmnet', 3, 4))),
         sites.num = ifelse(sites == 'Allsites', 1,
                            ifelse(sites == 'RFsites', 2, 3)),
         minCR.num = as.numeric(substr(minCR, 6, 6)),
         weight.num = ifelse(weight == 'none', 1, 
                             ifelse(weight == 'CR', 2, 3)),
         site.cr.num = substr(site.select.cr, 3, 3),
         site.cr.num = ifelse(site.cr.num == 'n', 0, site.cr.num),
         site.cr.num = as.numeric(site.cr.num),
         age.trans.num = ifelse( age.transform == 'none', 1, 2)
         ) |> 
  select(c(method.num, sites.num, minCR.num, weight.num, site.cr.num, age.trans.num, MAE))

fit <- rfPermute(
  formula = as.formula('MAE ~ .'), 
  data = temp,
  ntree = 100,
  replace = FALSE,
  num.rep = nrep,
  num.cores = 100
)

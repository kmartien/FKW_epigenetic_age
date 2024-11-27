library(ggplot2)
library(tidyverse)
load('data/age_and_methylation_data.rdata')

# CRC age distributions
p <- age.df %>%
  mutate(
    id = factor(swfsc.id, levels = .$swfsc.id[order(-age.best, -age.min, -age.max)]),
    Confidence = factor(age.confidence)
  ) |> 
  ggplot(aes(y = id)) +
  geom_segment(
    aes(x = age.min, xend = age.max, yend = id, color = Confidence), 
    linewidth = 2, alpha = 0.8
  ) +
  geom_point(
    aes(x = age.best, color = Confidence), 
    shape = 21, size = 3.5, fill = "white"
  ) + 
  scale_color_manual(values = conf.colors) +
  labs(x = "Age", y = NULL) + 
  scale_x_continuous(breaks = seq(0, max(age.df$age.max), 5)) +
  theme(
    text = element_text(size = 25),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.85, 0.85)
  )
jpeg(file = 'R/summaries/forPerth/CRC.age.distributions.jpg', width = 1000, height = 1000)
print(p)
dev.off()

conf.3plus.colors <- conf.colors
conf.3plus.colors[4] <- 'gray95'
p.CR3.plus <- age.df %>%
  mutate(
    id = factor(swfsc.id, levels = .$swfsc.id[order(-age.best, -age.min, -age.max)]),
    Confidence = factor(age.confidence)
  ) |> 
  ggplot(aes(y = id)) +
  geom_segment(
    aes(x = age.min, xend = age.max, yend = id, color = Confidence), 
    linewidth = 2, alpha = 0.8
  ) +
  geom_point(
    aes(x = age.best, color = Confidence), 
    shape = 21, size = 3.5, fill = "white"
  ) + 
  scale_color_manual(values = conf.3plus.colors) +
  labs(x = "Age", y = NULL) + 
  scale_x_continuous(breaks = seq(0, max(age.df$age.max), 5)) +
  theme(
    text = element_text(size = 25),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.85, 0.85)
  )
jpeg(file = 'R/summaries/forPerth/CRC.age.distributions.CR3plus.jpg', width = 1000, height = 1000)
print(p.CR3.plus)
dev.off()

conf.3plus.colors[3] <- 'gray95'
p.CR4.plus <- age.df %>%
  mutate(
    id = factor(swfsc.id, levels = .$swfsc.id[order(-age.best, -age.min, -age.max)]),
    Confidence = factor(age.confidence)
  ) |> 
  ggplot(aes(y = id)) +
  geom_segment(
    aes(x = age.min, xend = age.max, yend = id, color = Confidence), 
    linewidth = 2, alpha = 0.8
  ) +
  geom_point(
    aes(x = age.best, color = Confidence), 
    shape = 21, size = 3.5, fill = "white"
  ) + 
  scale_color_manual(values = conf.3plus.colors) +
  labs(x = "Age", y = NULL) + 
  scale_x_continuous(breaks = seq(0, max(age.df$age.max), 5)) +
  theme(
    text = element_text(size = 25),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.85, 0.85)
  )
jpeg(file = 'R/summaries/forPerth/CRC.age.distributions.CR4plus.jpg', width = 1500, height = 1500)
print(p.CR4.plus)
dev.off()

for(i in c(1:2)){
  jpeg(filename = paste0('R/summaries/forPerth/', names(base.plot)[[i]], '_base.scatterplot.plot.jpg'), width = 1200, height = 800, quality = 100)
  print(base.plot[[i]])
  dev.off()
}

jpeg(filename = paste0('R/summaries/forPerth/CR2.resample.results.jpg'), width = 1200, height = 1200)
age.dist[[2]]
dev.off()

age.dist.z0033906 <- 
  ggplot(filter(age.dist.gam.ranAgeMeth, swfsc.id == 'z0033906')) +
  geom_histogram(aes(x = age.pred, after_stat(density), fill = factor(age.confidence)), show.legend = FALSE) + 
  geom_vline(aes(xintercept = age.best, linewidth = 3), linetype = 'dashed', color = 'gray50') +
  geom_vline(aes(xintercept = lower, linewidth = 3)) +
  scale_fill_manual(values = conf.colors) +
  theme(text = element_text(size = 30)) +
  labs(x = 'Agepredicted', y = 'Density')
jpeg(file = 'R/summaries/forPerth/predicted.age.distribution.z0033906.jpg', height = 600, width = 1500)
age.dist.z0033906
dev.off()

jpeg(file = 'R/summaries/forPerth/predicted.age.distribution.wHDI.z0033906.jpg', height = 600, width = 1500)
age.dist.z0033906 + 
  ylim(c(0,0.1)) +
  geom_vline(
    aes(xintercept = lower, linewidth = 3),
    data =
      filter(age.dist.gam.ranAgeMeth, swfsc.id == 'z0033906') |>
        summarise(as.data.frame(rbind(HDInterval::hdi(age.pred))))
    ) +
  geom_vline(
    aes(xintercept = upper, linewidth = 3),
    data =
      filter(age.dist.gam.ranAgeMeth, swfsc.id == 'z0033906') |>
        summarise(as.data.frame(rbind(HDInterval::hdi(age.pred))))
    ) +
  geom_vline(
    aes(xintercept = age.best, linewidth = 3), 
    data = filter(age.dist.gam.ranAgeMeth, swfsc.id == 'z0033906'), 
    linetype = 'dashed', color = 'gray50')
dev.off()

# age.dist.gam.ranAgeMeth <- 
#   filter(ran.age.distributions, 
#          method == 'gam' & minCR == 'minCR4' & resample == 'ranAgeMeth') |> 
#   left_join(select(age.df, c(swfsc.id, age.confidence)))
# composite.age.dist <- lapply(1:nrow(age.df), function(i) {
#   crcAgeDist(age.df$swfsc.id[i], age.df, type = 'density') |> 
#     mutate(density = ifelse(age < age.df$age.min[i], 0, density),
#            density = ifelse(age > age.df$age.max[i], 0, density)) |> 
#     filter(age > age.df$age.min[i]-0.1 & age < age.df$age.max[i]+0.1)
# }) |>
#   bind_rows() |> 
#   left_join(select(age.df, c(swfsc.id, age.confidence)))

pred.GAM.ranAgeMeth <- left_join(
  pred.GAM.ranAgeMeth,
  age.dist.gam.ranAgeMeth |>  
  group_by(swfsc.id) |> 
  summarise(as.data.frame(rbind(HDInterval::hdi(age.pred)))) |> 
  rename(hdi.lower = lower, hdi.upper = upper),
  by = 'swfsc.id'
)

give.me.my.fucking.plots <- lapply(c(2,5), function(cr){ 
#  jpeg(file = paste0('R/summaries/forPerth/predicted.age.distributions.gam.ranAgeMeth.CR',cr,'.jpg'), height = 2000, width = 2000)
  age.dist <- ggplot(filter(age.dist.gam.ranAgeMeth, age.confidence == cr)) +
    geom_histogram(aes(x = age.pred, after_stat(density), fill = factor(age.confidence)), show.legend = FALSE) + 
    geom_vline(
      aes(xintercept = hdi.lower), 
      data = 
        filter(pred.GAM.ranAgeMeth, age.confidence == cr) , linewidth = 2
    ) +
    geom_vline(
      aes(xintercept = hdi.upper), 
      data = 
        filter(pred.GAM.ranAgeMeth, age.confidence == cr), linewidth = 2 
    ) +
    geom_line(aes(x = age, y = density), data = filter(composite.age.dist, age.confidence == cr), linewidth = 2, color = 'gray50') +
    geom_vline(
      aes(xintercept = age.best), 
      data = filter(pred.GAM.ranAgeMeth, age.confidence == cr), linetype = 'dashed', linewidth = 2, color = 'gray50') +
  scale_fill_manual(values = conf.colors) +
    labs(x = 'Predicted Age', y = 'Density') +
    theme(text = element_text(size = 25),
          legend.position = 'none') +
    facet_wrap(~swfsc.id, ncol = 5, scales = 'free')
    return(age.dist)
#  age.dist
#  dev.off()
})
jpeg(file = paste0('R/summaries/forPerth/predicted.age.distributions.gam.ranAgeMeth.CR2.jpg'), height = 500, width = 2000)
give.me.my.fucking.plots[[1]]
dev.off()
jpeg(file = paste0('R/summaries/forPerth/predicted.age.distributions.gam.ranAgeMeth.CR5.jpg'), height = 1000, width = 2000)
give.me.my.fucking.plots[[2]]
dev.off()

jpeg(file = 'R/summaries/forPerth/predicted.age.distributions.gam.ranAgeMeth.CR2.5.jpg', height = 2000, width = 2000)
do.call(grid.arrange, age.dist[[1]])
dev.off()

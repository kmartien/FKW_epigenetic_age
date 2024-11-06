library(ggplot2)
library(gridExtra)
library(viridis)

# plot loov results by min.CR
plot.loov.res <- function(loov.res, min.CR) {# expected argument is loov.res
  dat <- filter(loov.res, age.confidence >= min.CR)
  med.dev <- dat |> summarise(MAE = round(median(dev),2))
  cor.coeff <- round(cor.test(dat$age.best, dat$age.pred, method = "pearson")$estimate,2)
  loov.regress <- do.call(rbind,lapply(min.CR:5, function(cr){
    dat.filtered <- filter(dat, age.confidence == cr)
    regr <- lm(age.pred~age.best, data = dat.filtered)$coefficients
    names(regr) <- c("intercept", "slope")
    return(regr)
  }))
  fit.sum <- bind_cols(MAE = med.dev, Corr = cor.coeff, loov.regress)
  
  p.loov <- ggplot(dat, aes(x = age.best, y = age.pred)) +
    geom_point(aes(col = as.character(age.confidence)), size = 3) +
    # stat_smooth(aes(color = 'black', fill = NULL), 
    #             show.legend = FALSE,
    #             method = "lm", 
    #             formula = y ~ x, 
    #             geom = "smooth") +
    geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
    labs(x = "Age.best", y = "Predicted age") +
    annotation_custom(tableGrob(fit.sum[1, c('MAE', 'Corr')], theme = ttheme_minimal(base_size = 16), rows = NULL), xmin = 30, xmax = 40, ymin = 0, ymax = 10) +
    scale_color_manual(values = conf.colors, name = "Confidence") +
    xlim(0,40) + ylim(0,60) +
    theme_minimal() +
    theme(
      text = element_text(size = 20),
      legend.position = c(0.2, 0.8),
      plot.background = element_rect(color = "black", linewidth = 1)
    )
  # for(i in 1:nrow(fit.sum)){
  #   cr <- i+(min.CR-1)
  #   p.loov <- p.loov + 
  #     geom_abline(slope = fit.sum$slope[i], intercept = fit.sum$intercept[i], color = conf.palette[cr])
  # }
  return(list(p.loov = p.loov, fit.sum = fit.sum))
}

#scatterplot of residuals
plot.residuals <- function(loov.res, min.CR) {# expected argument is loov.res
  dat <- filter(loov.res, age.confidence >= min.CR)
  med.dev <- dat |> summarise(MAE = round(median(dev),2))
  cor.coeff <- round(cor.test(dat$age.best, dat$age.pred, method = "pearson")$estimate,2)
  # loov.regress <- do.call(rbind,lapply(min.CR:5, function(cr){
  #   dat.filtered <- filter(dat, age.confidence == cr)
  #   regr <- lm(age.pred~age.best, data = dat.filtered)$coefficients
  #   names(regr) <- c("intercept", "slope")
  #   return(regr)
  # }))
  fit.sum <- bind_cols(MAE = med.dev, Corr = cor.coeff)
  
  p.loov <- ggplot(dat, aes(x = age.best, y = resid)) +
    geom_point(aes(col = as.character(age.confidence)), size = 3) +
    stat_smooth(aes(color = 'black', fill = NULL), 
                show.legend = FALSE,
                method = "lm", 
                formula = y ~ x, 
                geom = "smooth") +
    geom_abline(slope = 0, color = "black", linewidth = 0.5, linetype = 2) +
    labs(x = "Age.best", y = "Residual") +
    annotation_custom(tableGrob(fit.sum[1, c('MAE', 'Corr')], theme = ttheme_minimal(base_size = 16), rows = NULL), xmin = 30, xmax = 40, ymin = 0, ymax = 10) +
    scale_color_manual(values = conf.colors, name = "Confidence") +
 #   xlim(0,40) + ylim(0,80) +
    theme_minimal() +
    theme(
      text = element_text(size = 20),
      legend.position = c(0.2, 0.2),
      plot.background = element_rect(color = "black", linewidth = 1)
    )
  return(p.loov)
}


rm(list = ls())
library(tidyverse)
library(sn)
library(swfscMisc)
library(parallel)
library(abind)
library(gridExtra)
source("R/misc_funcs.R")
source("R/data.prep/summarize_reads.R")
source("R/data.prep/calc_methylation.R")
set.seed(1)


# Load and format age data ------------------------------------------------

crc.data <- read.csv(
  'data_raw/Pseudorca_AgeEstimates_Simplified_2024MARv2.csv',
  na.strings = c(NA, 'NA', '')
  ) |> 
    select(c(ID.., Sex..G., LABID,Biopsy.., Date.of.biopsy, Best.age.point.estimate.at.time.of.sampling..estimated.in.2023,
             Best.age.point.estimate.confidence.rating..CR...1.low..5.high...estimated.in.2023, Minimum.age.point.at.time.of.sampling..estimated.in.2023,
             ALL.IDs..Maximum.age.estimate.2023, Pair..))
names(crc.data) <- c('crc.id', 'sex', 'swfsc.id', 'biopsy.id', 'date.biopsy', 'age.best', 'age.confidence', 'age.min', 'age.max', 'pair.id')

# split out and replicate ids
crc.data <- do.call(
  rbind,
  lapply(1:nrow(crc.data), function(i) {
    id <- strsplit(crc.data$swfsc.id[i], ",")[[1]] |> 
      stringi::stri_trim()
    df <- crc.data[rep(i, length(id)), ]
    df$swfsc.id <- id
    df
  })
)

# format labid
crc.data <- crc.data |> 
  mutate(
    swfsc.id = zero.pad(as.numeric(swfsc.id)),
    swfsc.id = ifelse(is.na(swfsc.id), NA, paste0("z0", swfsc.id))
  )


# format age class
age.class.map <- data.frame(
  class = c("calf", "juvenile", "sub.adult", "adult.f.sm", "adult.pm"),
  f.ages = c(2, 5, 9, max(crc.data$age.best), NA),
  m.ages = c(2, 5, 9, 24, max(crc.data$age.best))
)
crc.data <- crc.data |>  
  mutate(
    age.class = ifelse(
      sex == "Female", 
      as.numeric(cut(age.best, c(0, age.class.map$f), include = T, right = T)),
      as.numeric(cut(age.best, c(0, age.class.map$m), include = T, right = T))
    ),
    age.class = factor(
      age.class.map$class[age.class], 
      levels = age.class.map$class
    )
  )

crc.cols <- c(
  "crc.id", "swfsc.id", "biopsy.id", "date.biopsy", "sex",
  "age.best", "age.confidence", "age.min", "age.max", "age.class", "pair.id"
)
conf.map <- read_csv(
  'data_raw/CRC_confidence_rating_mappings.csv',
  name_repair = 'minimal',
  show_col_types = FALSE
) |> 
  column_to_rownames('person') |> 
  as.matrix()
age.df <- crc.data[, crc.cols] |> 
  mutate( 
    confidence.wt = unname(conf.map["sm", age.confidence]),
    date.biopsy = as.POSIXct(date.biopsy, format = "%d-%b-%y"),
    sex.num = as.numeric(factor(sex))
  )  |> 
  filter(!is.na(age.min)) |> 
  arrange(swfsc.id)
attributes(age.df$date.biopsy)$tzone <- "HST"


# change minimum age of z0178364 to 10 since KMR GLG age was 12.5 and CRC best
# age is 10
age.df$age.min[age.df$swfsc.id == 'z0178364'] <- 10


# fit Skew Normal parameters
age.sn.optim <- mclapply(1:nrow(age.df), function(i) {
  sn.params(
    age.df$age.best[i],
    age.df$age.min[i],
    age.df$age.max[i],
    p = 0.975,
    shape.const = 4
  )
}, mc.cores = 10)
age.sn.params <- t(sapply(age.sn.optim, function(x) x$par))
colnames(age.sn.params) <- c("sn.location", "sn.scale", "sn.shape")
age.df <- cbind(age.df, age.sn.params) |> 
  mutate(
    age.range = age.max - age.min,
    age.var = sn.scale ^ 2 * (1 - ((2 * swfscMisc::sn.delta(sn.shape) ^ 2) / pi)),
    dens.unif = 1 / age.range,
    wtd.unif = (1 - confidence.wt) / age.range,
    sn.age.offset = age.best - (sn.scale * swfscMisc::sn.m0(sn.shape))
  ) 

age.df <- arrange(age.df, swfsc.id)


# CRC confidence colors ---------------------------------------------------

conf.colors <- setNames(
  colorRampPalette(viridis::magma(30, direction = -1)[-(1:5)])(5),
  5:1
)
conf.colors[5] <- 'gray40'


# Load and format methylation data ----------------------------------------

meth.dat <- summarize.reads(min.cov = 100)

mean.cov <- do.call("cbind",lapply(meth.dat,function(s){s$CpG.sum$coverage$avg}))
rownames(mean.cov) <- meth.dat$TET2$CpG.sum$coverage$id
colnames(mean.cov) <- names(meth.dat)
write.csv(mean.cov, "results_raw/mean.coverage.Aligned.to.targets.with.primers.csv")

# Calculate conversion efficiency -----------------------------------------------
pct.conversion <- left_join(
  tibble(id = meth.dat[[1]]$non.CpG.sum$freq.meth$id,
         tot.cov = colSums(do.call('rbind', lapply(meth.dat, function(s){
           last.pos <- dim(s$non.CpG.sum$freq.meth)[2]
           if (!is.null(last.pos)){
             tot.cov <- rowSums(s$non.CpG.sum$coverage[,2:last.pos],na.rm = TRUE) - rowSums(s$non.CpG.sum$errors[,2:last.pos],na.rm = TRUE)
           }
         })))),
  
  tibble(id = meth.dat[[1]]$non.CpG.sum$freq.meth$id,
         tot.Ts = colSums(do.call('rbind',lapply(meth.dat, function(s){
           last.pos <- dim(s$non.CpG.sum$freq.meth)[2]
           if (!is.null(last.pos)){
             tot.Ts <- (rowSums(s$non.CpG.sum$coverage[,2:last.pos],na.rm = TRUE)-rowSums(s$non.CpG.sum$freq.meth[,2:last.pos],na.rm = TRUE)-rowSums(s$non.CpG.sum$errors[,2:last.pos],na.rm = TRUE))
           }
         })))),
  by = 'id'
) |> 
  mutate(conversion = tot.Ts / tot.cov)

save(pct.conversion, file = "data/pct.conversion.rda")

# Logit transformed methylation -------------------------------------------------------

logit.meth <- do.call('cbind',lapply(1:length(meth.dat), function(amp){
  s <- meth.dat[[amp]]
  last.pos <- dim(s$CpG.sum$freq.meth)[2]
  meth <- s$CpG.sum$freq.meth[,2:last.pos]
  corrected.cov <- s$CpG.sum$coverage[,2:last.pos] - s$CpG.sum$errors[,2:last.pos]
  x <- meth/(corrected.cov)
  x <- 1-((1-x)/pct.conversion$conversion)
  rownames(x) <- s$CpG.sum$freq.meth[,1]
  colnames(x) <- paste(names(meth.dat)[amp],zero.pad(as.numeric(colnames(meth))),sep="_")
  for (j in 1:ncol(x)){
    col.min.nonzero <- min(x[which(x[,j] > 0),j])
    x[which(x[,j] <= 0),j] <- col.min.nonzero * 0.5
  }
  return(log(x/(1-x)))
}))

# IDs with good coverage and no missing sites -----------------------------

cpg <- do.call(
  rbind,
  lapply(names(meth.dat), function(locus) {
    do.call(
      rbind,
      lapply(names(meth.dat[[locus]]), function(type) {
        x <- meth.dat[[locus]][[type]]
        if(is.null(x)) return(NULL)
        x$coverage |> 
          select(-avg) |> 
          pivot_longer(-id, names_to = "site", values_to = "coverage") |> 
          left_join(
            pivot_longer(x$freq.meth, -id, names_to = "site", values_to = "freq.meth"),
            by = c("id", "site")
          ) |> 
          left_join(
            pivot_longer(x$errors, -id, names_to = "site", values_to = "errors"),
            by = c("id", "site")
          ) |>
          mutate(locus = locus, type = type)
      })
    )
  })
) |>   
  mutate(
    site = as.numeric(site),
    loc.site = paste0(locus, "_", zero.pad(site)),
    id.site = paste0(id, "_", loc.site),
    corrected.cov = coverage - errors,
    pct.meth = freq.meth / corrected.cov,
    type = gsub(".sum", "", type)
  ) |> 
  arrange(type, id, locus, site) |> 
  filter(type == 'CpG') |> 
  left_join(
    select(pct.conversion, c(id, conversion)),
    by = 'id'
    )

ids.to.keep <- cpg |> 
  filter(id %in% age.df$swfsc.id & type == "CpG") |> 
  group_by(id) |> 
  summarize(
    is.complete = !any(is.na(freq.meth)),
    median.cov = median(coverage, na.rm = TRUE), 
    .groups = "drop"
  ) |> 
  # eliminate low coverage samples and z0091086 because it ends up with coverage < 100 at one of the CpG site we keep
  filter(median.cov >= 1000 & id != 'z0091086') |> 
  pull("id")


# Select sites ----------------------------------------------

sites.to.keep <- cpg |> 
  filter(id %in% ids.to.keep) |> 
  group_by(loc.site) |> 
  summarize(
    cov.median = median(corrected.cov, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  filter(cov.median >= 1000) |> 
  pull("loc.site")


# Plot age distributions --------------------------------------------------

graphics.off()
pdf(paste0("results_raw/summary - CRC age priors.pdf"), height = 10, width = 10)

p <- filter(age.df, swfsc.id %in% ids.to.keep) %>%
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
    text = element_text(size = 10),
    legend.position = "top"
  )
print(p)

p <- age.df |> 
  mutate(Confidence = factor(age.confidence)) |> 
  ggplot(aes(age.best, age.range)) +
  geom_point(aes(fill = Confidence), color = "white", shape = 21, size = 4) +
  scale_fill_manual(values = conf.colors) +
  labs(x = "CRC best age", y = "Maximum - minimum age") +
  theme(legend.position = "top")
print(p)


for(cr.df in split(age.df, age.df$age.confidence)) {
  conf <- unique(cr.df$age.confidence)
  p <- lapply(1:nrow(cr.df), function(i) {
    crcAgeDist(cr.df$swfsc.id[i], cr.df, type = 'density')
  }) |>
    bind_rows() |>
    ggplot() +
    geom_vline(aes(xintercept = age.min), linetype = 'dashed', data = cr.df) +
    geom_vline(aes(xintercept = age.max), linetype = 'dashed', data = cr.df) +
    geom_vline(aes(xintercept = age.best), data = cr.df) +
    geom_line(
      aes(x = age, y = density),
      color = conf.colors[as.character(conf)],
      linewidth = 1.5
    ) +
    labs(x = 'Age', y = 'Density', title = paste0("Confidence rating = ", conf)) +
    facet_wrap(~ swfsc.id, scales = 'free_y') +
    theme_minimal() +
    theme(
      axis.line.x = element_line(),
      axis.line.y = element_line()
    )
  print(p)
}

dev.off()

# Paired sample age difference --------------------------------------------

pair.df <- age.df |>
  filter(!is.na(age.df$pair.id) & swfsc.id %in% ids.to.keep)

num.samples <- 1000000
pair.df <- split(pair.df, pair.df$pair.id) |>
  lapply(function(df) {
    if(n_distinct(df$swfsc.id) < 2) NULL else {
      combn(df$swfsc.id, 2, function(id) {
        pairs <- data.frame(swfsc.id = id, pair.id = unique(df$pair.id)) |>
          left_join(select(age.df, swfsc.id, date.biopsy), by = 'swfsc.id') |>
          arrange(desc(date.biopsy))
        
        older <- pairs$swfsc.id[1]
        younger <- pairs$swfsc.id[2]
        
        older.sn.dp <- age.df |>
          filter(swfsc.id == older) |>
          select(sn.location, sn.scale, sn.shape) |>
          unlist()
        younger.sn.dp <- age.df |>
          filter(swfsc.id == younger) |>
          select(sn.location, sn.scale, sn.shape) |>
          unlist()
        pair.diff <- sn::rsn(num.samples, dp = older.sn.dp) - sn::rsn(num.samples, dp = younger.sn.dp)
        
        data.frame(
          older = older,
          younger = younger,
          pair.id = unique(df$pair.id),
          year.diff = as.numeric(difftime(pairs$date.biopsy[1], pairs$date.biopsy[2], 'days')) / 365.25,
          mean = mean(pair.diff),
          sd = sd(pair.diff)
        )
      }, simplify = FALSE) |>
        bind_rows()
    }
  }) |>
  bind_rows()

# Save data ---------------------------------------------------------------

save(
  age.df, conf.colors, cpg, meth.dat,
  ids.to.keep, sites.to.keep, pair.df, logit.meth,
  file = "data/age_and_methylation_data.rdata"
)



df <- logit.meth.normal.params |> 
  mutate(
    meth.lci = qnorm(0.025, mean.logit, sd.logit),
    meth.uci = qnorm(0.975, mean.logit, sd.logit)
  ) |> 
  left_join(
    select(age.df, age.best, age.min, age.max, age.confidence, swfsc.id),
    by = 'swfsc.id'
  ) |> 
  mutate(age.confidence = factor(age.confidence))

graphics.off()
pdf('summary - logit meth by age.pdf')
for(i in sort(sites.to.keep)) {
  g <- df |> 
    filter(loc.site == i) |> 
    ggplot() +
    geom_vline(
      aes(xintercept = median(mean.logit)), 
      color = 'gray10', alpha = 0.5, linetype = 'dashed'
    ) +
    geom_hline(
      aes(yintercept = median(age.df$age.best)),
      color = 'gray10', alpha = 0.5, linetype = 'dashed'
    ) +
    geom_segment(
      aes(x = meth.lci, xend = meth.uci, y = age.best, yend = age.best, color = age.confidence),
      alpha = 0.6, linewidth = 0.2
    ) +
    geom_segment(
      aes(x = mean.logit, xend = mean.logit, y = age.min, yend = age.max, color = age.confidence),
      alpha = 0.6, linewidth = 0.2
    ) +
    geom_point(
      aes(mean.logit, age.best, fill = age.confidence), 
      color = 'white', shape = 21, size = 3
    ) +
    scale_fill_manual(values = conf.colors) +
    scale_color_manual(values = conf.colors) +
    labs(x = 'logit(Pr(meth))', y = 'CRC age', title = i) + 
    theme(
      legend.position = 'top',
      legend.title = element_blank()
    )
  print(g)
}
dev.off()

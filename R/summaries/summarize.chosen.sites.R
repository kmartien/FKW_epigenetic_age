all.chosen.sites <- bind_rows(
  readRDS('R/glmnet_tuning/glmnet.chosen.sites.rds') |> 
    rename(loc.site = sites) |> 
    filter(count >= 500 & loc.site != '(Intercept)') |> 
    mutate(params = paste('glmnet', cr, wt, age.transform, sep = '_'),
           chosen = 'X') |> 
    select(c(loc.site, params, chosen)),
  readRDS('R/rf_tuning/rf_chosen_sites.rds') |> 
    filter(pval <= 0.05) |> 
    mutate(params = paste('rf', cr, wt, age.transform, sep = '_'),
           chosen = 'X') |> 
    select(c(loc.site, params, chosen))
) |> 
  pivot_wider(names_from = params, values_from = chosen) 

write.csv(all.chosen.sites, file = 'results_raw/chosen.site.summary.csv')

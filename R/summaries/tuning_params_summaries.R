library(tidyverse)

fnames <- tibble(name = list.files('R/rf_tuning', pattern = 'optim')) |> 
  mutate(len = nchar(name)) |> 
  filter(len < 30) |> 
  mutate(model = substr(name, 22, nchar(name) - 4))

rf_param_sum <- lapply(1:nrow(fnames), function(i){
  load(paste0('R/rf_tuning/', fnames$name[i]))
  return(tibble(model = fnames$model[i], rf.params))
}) |> bind_rows()
  

fnames <- tibble(name = list.files('R/glmnet_tuning', pattern = 'optim')) |> 
  mutate(model = substr(name, 18, nchar(name) - 4))

glmnet_param_sum <- lapply(1:nrow(fnames), function(i){
  optim_params <- readRDS(paste0('R/glmnet_tuning/', fnames$name[i]))
  return(tibble(model = fnames$model[i], alpha = optim_params$par))
}) |> bind_rows()

fnames <- tibble(name = list.files('R/svm', pattern = 'tuning')) |> 
  mutate(model = substr(name, 17, nchar(name) - 4))

svm_param_sum <- lapply(1:nrow(fnames), function(i){
  load(paste0('R/svm/', fnames$name[i]))
  return(tibble(model = fnames$model[i], tune.obj$best.parameters))
}) |> bind_rows()

## simple environmental sampler function. 
## Now this is re-calculated for each year (so last years temp doesn't line up with 
## the lagged temp in current year)
sampling_env <- function(iteration, env_params, start_year = 2023) {
  
  lags_sim <- matrix(c(1:(env_params$lags + 1)), nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  
  ### get values for window function
  # mon <- seq(from = 0, to = 0.95, length.out = 12)
  # end <- (start_year + iteration) + mon[7]
  # start <- end - (1/12 * (env_params$lags + 1))
  # 
  # pet_sim <- rev(as.numeric(window(env_params[[grep("pet", names(env_params), value = T)]], start, end))) %>%
  #   matrix(., nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  # pr_sim <- rev(as.numeric(window(env_params[[grep("pr", names(env_params), value = T)]], start, end))) %>%
  #   matrix(., nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  # tas_sim <- rev(as.numeric(window(env_params[[grep("tas", names(env_params), value = T)]], start, end))) %>%
  #   matrix(., nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  
  
  return(list(lags = lags_sim,
              precip = env_params$pr,
              pet = env_params$pet,
              shading = env_params$shading
  ))
}



FLM_clim_predict <- function(model, lag_vec, 
                             precip_vec, pet_vec, 
                             shading) {
  
  ## dummy datalist, won't be using the predictions for n_stems : year_t0 but these need to be provided for predict function
  new_data <- list(
    ln_stems_t0 = 2,  
    population = "CR",
    tot_shading_t0 = 0,
    soil_depth = 0,
    year_t0 = 2016,
    tot_shading_m = matrix(rep(shading, length(lag_vec)), nrow = 1),
    lags = matrix(lag_vec, nrow = 1),
    # tas_scaledcovar = matrix(temp_vec, nrow = 1),
    pr_scaledcovar = matrix(precip_vec, nrow = 1),
    pet_scaledcovar = matrix(pet_vec, nrow = 1)
  )
  
  pt <- mgcv::predict.gam(model, new_data, type = "terms")
  
  
  
  return(
    sum(pt[, grep("scaledcovar", attributes(pt)[[2]][[2]], value = T)])
  )
}

run_ipm <- function(params, env_params, 
                    n_it = n_it, U, L, n){
  
  init_ipm("general", "di", "stoch", "param") %>%
    define_kernel(
      name      = "P",
      family    = "CC",
      formula   = s * g * d_stems,
      
      s         =  plogis(s_linear),
      s_linear  =  s_int + s_stems * stems_1 + s_site + 
        FLM_clim_predict(model = surv_mod, ### spline model prediction
                         lag_vec = lags,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading
        ),
      g         =  dnorm(stems_2, g_mu, grow_sd),
      g_mu      =  g_int + g_stems * stems_1 + g_site +
        FLM_clim_predict(model = grow_mod,
                         lag_vec = lags,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading), ### model prediction
      
      data_list     = params,
      states        = list(c("stems")),
      
      uses_par_sets = F,
      # par_set_indices = list(loc = loc),
      
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", "g")
    ) %>%
    define_kernel(
      name = "F_to_SDL",
      family = "CD",
      formula = fp * pabort * n_seeds * surv1_seed * germ * d_stems,
      
      fp          = plogis(fp_linear),
      fp_linear   = fp_int + fp_stems * stems_1 + fp_site +
        FLM_clim_predict(model = pflower_mod,
                         lag_vec = lags,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      pabort     = plogis(ab_linear),
      ab_linear   = ab_int + ab_stems * stems_1 + ab_site +
        FLM_clim_predict(model = pabort_mod,
                         lag_vec = lags,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      n_seeds = ifelse(n_seeds_mod > 515, 515, n_seeds_mod),
      n_seeds_mod     = exp(n_seeds_linear),
      n_seeds_linear = ns_int + ns_stems * stems_1 + ns_site +
        FLM_clim_predict(model = nseed_mod,
                         lag_vec = lags,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      germ        = germ_mean, 
      surv1_seed  = seed_surv1,
      
      data_list = params,
      states = list(c("stems", "sdl")),
      
      uses_par_sets = F,
      # par_set_indices = list(loc = loc),
      
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name = "F_to_SB1",
      family = "CD",
      formula = fp * pabort * n_seeds * surv1_seed * (1-germ) * d_stems,
      
      fp          = plogis(fp_linear),
      fp_linear   = fp_int + fp_stems * stems_1 + fp_site +
        FLM_clim_predict(model = pflower_mod,
                         lag_vec = lags,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      pabort     = plogis(ab_linear),
      ab_linear   = ab_int + ab_stems * stems_1 + ab_site +
        FLM_clim_predict(model = pabort_mod,
                         lag_vec = lags,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      n_seeds = ifelse(n_seeds_mod > 645, 645, n_seeds_mod),
      n_seeds_mod     = exp(n_seeds_linear),
      n_seeds_linear = ns_int + ns_stems * stems_1 + ns_site +
        FLM_clim_predict(model = nseed_mod,
                         lag_vec = lags,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      germ        = germ_mean, 
      surv1_seed  = seed_surv1,
      
      data_list = params,
      states = list(c("stems", "sb1")),
      
      uses_par_sets = F,
      # par_set_indices = list(loc = loc),
      
      evict_cor = FALSE
      
    ) %>%
    define_kernel(
      name = "SB1_to_SB2",
      family = "DD",
      formula = (1-germ) * surv2_seed,
      
      germ = germ_mean,
      surv2_seed = seed_surv2,
      
      data_list = params,
      states = list(c("sb1", "sb2")),
      
      uses_par_sets = FALSE,
      
      evict_cor = FALSE
      
    ) %>%
    define_kernel(
      name = "SB1_to_SDL",
      family = "DD",
      formula = germ * surv2_seed,
      
      germ = germ_mean,
      surv2_seed = seed_surv2,
      
      data_list = params,
      states = list(c("sb1", "sdl")),
      
      uses_par_sets = FALSE,
      
      evict_cor = FALSE
      
    ) %>%
    define_kernel(
      name = "SB2_to_SDL",
      family = "DD",
      formula = germ * surv3_seed,
      
      germ       = germ_mean,
      surv3_seed = seed_surv3,
      
      data_list = params,
      states = list(c("sb2", "sdl")),
      
      uses_par_sets = FALSE,
      
      evict_cor = FALSE
      
    ) %>%
    define_kernel(
      name = "SDL_to_Plant",
      family = "DC",
      formula = sdl_surv * d_size * d_stems,
      
      sdl_surv     = plogis(sdl_s_linear),
      sdl_s_linear = sdl_s_int,
      d_size       = dnorm(stems_2, sdl_d_linear, sdl_size_d_sd),
      sdl_d_linear = sdl_d_int,
      
      data_list = params,
      states = list(c("sdl", "stems")),
      
      uses_par_sets = F,
      # par_set_indices = list(loc = loc),
      
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "d_size")
      
    ) %>%
    define_impl(
      list(
        "P" =        list(int_rule = 'midpoint',
                            state_start = "stems",
                            state_end = "stems"),
        "F_to_SDL" =     list(int_rule = 'midpoint',
                                state_start = "stems",
                                state_end = "sdl"),
        "F_to_SB1" =     list(int_rule = 'midpoint',
                                state_start = "stems",
                                state_end = "sb1"),
        "SB1_to_SB2" =   list(int_rule = "midpoint",
                            state_start = "sb1",
                            state_end = "sb2"),
        "SB1_to_SDL" =   list(int_rule = "midpoint",
                            state_start = "sb1",
                            state_end = "sdl"),
        "SB2_to_SDL" =   list(int_rule = "midpoint",
                            state_start = "sb2",
                            state_end = "sdl"),
        "SDL_to_Plant" = list(int_rule = "midpoint",
                                state_start = "sdl",
                                state_end = "stems")
      )
    ) %>%
    define_domains(
      stems = c(L, U, n)
    ) %>%
    define_env_state(
      env_covs = sample_env(iteration = t, env_params = env_params),
      data_list = list(env_params = env_params,
                       sample_env = sampling_env)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_stems = runif(n),
        n_sb1 = 20,
        n_sb2 = 20,
        n_sdl = 20
      )
    ) %>%
    make_ipm(
      iterate = T,
      iterations = n_it,
      # kernel_seq = rep(locality, n_it),
      usr_funs = list(FLM_clim_predict = FLM_clim_predict),
      return_sub_kernels = T
    )
  
}

ipm_loop <- function(yr, loc, shading, params,
                       climate_df,
                     n_it, U = U, L = L, n = n) {
  #print(i)
  
  # environmental params
  env_params <- list(
    pet = climate_df$pet_scaledcovar[which(climate_df$year_t0 == yr & climate_df$localities == loc),],
    pr = climate_df$pr_scaledcovar[which(climate_df$year_t0 == yr & climate_df$localities == loc),],
    lags = 24,
    shading = shading
    )
  
  # This line changes the name of the site specific parameters into general names to use in the IPM
  names(params) <- gsub(pattern = paste0("_", loc), replacement = "", names(params))
    
  ipm <- run_ipm(params = params, env_params = env_params, 
                 n_it = n_it, U = U, L = L, n = n)
  
  df1 <- data.frame(year_t0 = yr,
                    locality = loc,
                    shading = shading,
                    lambda = ipmr::lambda(ipm))
  
  return(df1)
}



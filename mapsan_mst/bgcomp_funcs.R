### Functions to Perform Bayesian g-computation for multiple binary outcomes
### David Holcomb
### 2024-12-30


### Instructions
# 'source' this script to load these functions into your R session:
#         source("bgcomp_funcs.R")
#
# These functions require the following packages to be loaded:
#         ## Stan Packages
#         library(brms)
#         library(cmdstanr)
#         library(tidybayes)
# 
#         ## Basic Packages
#         library(tidyverse)
# 
# `cmdstan` and `cmdstanr` should be installed on your computer before installing `brms`
# installing `cmdstan` requires a properly configured C++ toolchain:
#       https://mc-stan.org/docs/cmdstan-guide/installation.html#cpp-toolchain
#
# then install`cmdstanr` and `cmdstan` following these instructions:
#       https://mc-stan.org/cmdstanr/articles/cmdstanr.html
#       



### helper functions
## determine if variable is binary:
is.binary <- function(x) { all(na.omit(x) %in% 0:1) }

## define NOT IN
`%notin%` = function(x, y) { !(x %in% y) }



### BINARY EXPOSURE posterior predictions
bgcomp.bin <- function(fit, data,
                       exposure, outcome, group){
  # row-wise posterior predictive expectations from fitted brms model object
  
  if(is.null(group)){
    
    df_pred <- data
    
  } else {
    
    ## add pooled outcome data (new group)
    df_pred <- data %>%
      mutate(!!group := "pooled",
             !!outcome := NA) %>%
      distinct() %>%
      bind_rows(data)
    
  }
  
  
  ## prediction data frame
  dat_pred <- bind_rows(df_pred %>%
                          mutate(expose_obs = .data[[exposure]],
                                 {{exposure}} := expose_obs,
                                 contrast = "p_observed"),
                        df_pred %>%
                          mutate(expose_obs = .data[[exposure]],
                                 {{exposure}} := 0,
                                 contrast = "p_unexposed"),
                        df_pred %>%
                          mutate(expose_obs = .data[[exposure]],
                                 {{exposure}} := 1,
                                 contrast = "p_exposed")
  )
  
  
  ## posterior predictive expectation
  ppe <- epred_draws(object = fit,
                     newdata = dat_pred,
                     allow_new_levels = TRUE,
                     sample_new_levels = "uncertainty") %>%
    ungroup() %>%
    select(-c(.row, .chain, .iteration, {{exposure}})) %>%
    pivot_wider(names_from = contrast,
                values_from = .epred) %>%
    mutate(RD = p_exposed - p_unexposed,          # marginal risk difference
           RR = p_exposed / p_unexposed,          # marginal risk ratio
           OR = (p_exposed / (1 - p_exposed)) /   # marginal odds ratio
             (p_unexposed / (1 - p_unexposed)),
           support = "binary")
  
  
  return(ppe)
  
}



### CONTINUOUS EXPOSURE posterior predictions
bgcomp.cont <- function(fit, data,
                        exposure, outcome, group,
                        eps = 1e-4){
  
  # row-wise posterior predictive expectations from fitted brms model object
  
  if(is.null(group)){
    
    df_pred <- data
    
  } else {
    
    ## add pooled outcome data (new group)
    df_pred <- data %>%
      mutate(!!group := "pooled",
             !!outcome := NA) %>%
      distinct() %>%
      bind_rows(data)
    
  }
  
  ## prediction data frame
  dat_pred <- bind_rows(df_pred %>%
                          mutate(expose_obs = .data[[exposure]],
                                 {{exposure}} := expose_obs,
                                 contrast = "p_observed"),
                        df_pred %>%
                          mutate(expose_obs = .data[[exposure]],
                                 {{exposure}} := expose_obs - eps/2,
                                 contrast = "p_unexposed"),
                        df_pred %>%
                          mutate(expose_obs = .data[[exposure]],
                                 {{exposure}} := expose_obs + eps/2,
                                 contrast = "p_exposed")
  )
  
  
  ## posterior predictive expectation
  ppe <- epred_draws(object = fit,
                     newdata = dat_pred,
                     allow_new_levels = TRUE,
                     sample_new_levels = "uncertainty") %>%
    ungroup() %>%
    select(-c(.row, .chain, .iteration, {{exposure}})) %>%
    pivot_wider(names_from = contrast,
                values_from = .epred) %>%
    mutate(RD = (p_exposed - p_unexposed) / eps,  # marginal slope (dY/dX)--risk difference for 1 unit increase
           RR = (p_observed + RD) / p_observed,   # marginal risk ratio (1 unit increase)
           OR = ((p_observed + RD) /              # marginal odds ratio (1 unit increase)
                   (1 - (p_observed + RD))) / 
             (p_observed / (1 - p_observed)),
           support = "continuous")
  
  
  return(ppe)
  
}



### Overall function to fit model and perform g-computation
bgcomp <- function(data,                  # input data frame
                   formula,               # regression model formula (brms style)
                   summary = TRUE,        # summarize predictions across observations?
                   exposure,              # exposure variable name (string)
                   outcome = "infect",    # outcome variable name (string) 
                   group = "path",        # grouping variable name (string)
                   eps = 1e-4,            # epsilon size for finite difference (scalar)
                   prior = NULL,          # brms model priors (default brms improper flat priors)
                   sample_prior = "no",   # set to "yes" to sample prior predictive distribution (requires proper priors)
                   file = NULL,           # provide file path & name to save model object as RDS
                   file_refit = "always", # if model is saved, overwrite ("always, default), reuse ("never"), 
                                          #    or only overwrite if model changes ("on_change")
                   chains = 4,            # independent MCMC chains
                   iter = 2000,           # total MCMC iterations per chain (including warmup)
                   warmup = floor(iter/2),# warmup iterations to discard
                   cores = getOption("mc.cores", 1)){
  # data should already be filtered for outcome & complete cases
  # currently only accepts binary outcomes
  # must include a unique exposure id column if multiple exposures per subject
  
  
  ### fit the regression model
  
  ## QC formula
  form <- as.formula(formula)
  
  ## Determine exposure type
  expose_bin <- data %>%
    select({{exposure}}) %>%
    pull() %>%
    is.binary()
  
  if(expose_bin){
    data <- data %>%
      mutate({{exposure}} := as.integer(.data[[exposure]]))
  }
  
  ## brms bayesian hierarchical model
  fit <- brm(data = data,
             formula = form,
             family = bernoulli(),  # currently only set up for binary outcomes
             prior = prior,
             sample_prior = sample_prior,
             backend = "cmdstanr",
             file = file,
             file_refit = file_refit,
             chains = chains,
             iter = iter,
             warmup = warmup,
             cores = cores)
  
  
  ### posterior predictions
  
  ## Make posterior predictions
  
  if(expose_bin){
    # binary exposure
    ppe <- bgcomp.bin(fit = fit, data = data,
                      exposure = exposure, outcome = outcome, group = group)
    
  } else {
    # continuous exposure
    ppe <- bgcomp.cont(fit = fit, data = data,
                       exposure = exposure, outcome = outcome, group = group,
                       eps = eps)
  }
  
  
  
  ### Summarise and return 
  if(summary){
    
    if(is.null(group)){
      
      group <- "path"
      
      ppe <- ppe %>%
        mutate(path = "single")
      
    }
        
    ppe <- ppe %>%
      select({{group}},
             p_observed, p_unexposed, p_exposed,
             RD, RR, OR) %>%
      pivot_longer(c(p_observed, p_unexposed, p_exposed,
                     RD, RR, OR),
                   names_to = "param",
                   values_to = "pred") %>%
      group_by(across({{group}}), param) %>%
      summarise(avg = mean(pred),
                sd = sd(pred),
                iqr = IQR(pred),
                mad = mad(pred),
                p2.5 = quantile(pred, prob = 0.025),
                p5 = quantile(pred, prob = 0.05),
                p10 = quantile(pred, prob = 0.1),
                p25 = quantile(pred, prob = 0.25),
                p50 = quantile(pred, prob = 0.5),
                p75 = quantile(pred, prob = 0.75),
                p90 = quantile(pred, prob = 0.9),
                p95 = quantile(pred, prob = 0.95),
                p97.5 = quantile(pred, prob = 0.975)) %>%
      ungroup()
    
  }
  
  return(ppe)
  
}



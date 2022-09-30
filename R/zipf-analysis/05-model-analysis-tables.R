# This script prepares/analyzes data comparing correct/incorrect model specifications


### load libraries ---
library(tidyverse)

### load data ----
dtm_metrics <- read_rds("data-derived/zipf-analysis/dtm-metrics.rds")

model_metrics <- read_rds("data-derived/zipf-analysis/model-metrics.rds")

### format model metrics table ---
model_metrics <- 
  model_metrics |>
  mutate(
    wrong_k = Nk != est_k,
    wrong_eta = beta_sum != est_eta_sum & ! est_flat_eta,
    wrong_alpha = alpha_sum != est_alpha_sum,
    small_k = Nk > est_k,
    big_k = Nk < est_k,
    small_eta = beta_sum > est_eta_sum,
    big_eta = beta_sum < est_eta_sum,
    small_alpha = alpha_sum > est_alpha_sum,
    big_alpha = alpha_sum < est_alpha_sum
  ) |>
  mutate(
    correct_model = ! est_flat_eta & 
      ! wrong_k &
      ! wrong_eta &
      ! wrong_alpha,
    id = (1:(nrow(model_metrics) / 5)) |>
      rep(5) |>
      sort()
  )

### Create tables of pairs based on correct/incorrect spec by hyperparam ----

# group by model ID
matched_model_metrics <- 
  model_metrics |>
  by(INDICES = model_metrics$id, function(x) x)


# If K is too big
nk_big <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_k != x$Nk
      
      if (x$est_k[x$correct_model] < x$est_k[indicator]) {
        out <- tibble(
          id = x$id[1],
          ll_correct = x$ll[x$correct_model],
          ll_est = x$ll[indicator],
          r2_correct = x$r2[x$correct_model],
          r2_est = x$r2[indicator],
          mean_coherence_correct = x$mean_coherence[x$correct_model],
          mean_coherence_est = x$mean_coherence[indicator],
          var_coherence_correct = x$var_coherence[x$correct_model],
          var_coherence_est = x$var_coherence[indicator],
          mean_prevalence_correct = x$mean_prevalence[x$correct_model],
          mean_prevalence_est = x$mean_prevalence[indicator],
          var_prevalence_correct = x$var_prevalence[x$correct_model],
          var_prevalence_est = x$var_prevalence[indicator],
          geweke_correct = x$geweke_stat[x$correct_model],
          geweke_est = x$geweke_stat[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()

# if K is too small
nk_small <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_k != x$Nk
      
      if (x$est_k[x$correct_model] > x$est_k[indicator]) {
        out <- tibble(
          id = x$id[1],
          ll_correct = x$ll[x$correct_model],
          ll_est = x$ll[indicator],
          r2_correct = x$r2[x$correct_model],
          r2_est = x$r2[indicator],
          mean_coherence_correct = x$mean_coherence[x$correct_model],
          mean_coherence_est = x$mean_coherence[indicator],
          var_coherence_correct = x$var_coherence[x$correct_model],
          var_coherence_est = x$var_coherence[indicator],
          mean_prevalence_correct = x$mean_prevalence[x$correct_model],
          mean_prevalence_est = x$mean_prevalence[indicator],
          var_prevalence_correct = x$var_prevalence[x$correct_model],
          var_prevalence_est = x$var_prevalence[indicator],
          geweke_correct = x$geweke_stat[x$correct_model],
          geweke_est = x$geweke_stat[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()


# if eta_sum is too big
eta_sum_big <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_eta_sum != x$beta_sum
      
      if (x$est_eta_sum[x$correct_model] < x$est_eta_sum[indicator]) {
        out <- tibble(
          id = x$id[1],
          ll_correct = x$ll[x$correct_model],
          ll_est = x$ll[indicator],
          r2_correct = x$r2[x$correct_model],
          r2_est = x$r2[indicator],
          mean_coherence_correct = x$mean_coherence[x$correct_model],
          mean_coherence_est = x$mean_coherence[indicator],
          var_coherence_correct = x$var_coherence[x$correct_model],
          var_coherence_est = x$var_coherence[indicator],
          mean_prevalence_correct = x$mean_prevalence[x$correct_model],
          mean_prevalence_est = x$mean_prevalence[indicator],
          var_prevalence_correct = x$var_prevalence[x$correct_model],
          var_prevalence_est = x$var_prevalence[indicator],
          geweke_correct = x$geweke_stat[x$correct_model],
          geweke_est = x$geweke_stat[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()


# if eta_sum is too small
eta_sum_small <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_eta_sum != x$beta_sum
      
      if (x$est_eta_sum[x$correct_model] > x$est_eta_sum[indicator]) {
        out <- tibble(
          id = x$id[1],
          ll_correct = x$ll[x$correct_model],
          ll_est = x$ll[indicator],
          r2_correct = x$r2[x$correct_model],
          r2_est = x$r2[indicator],
          mean_coherence_correct = x$mean_coherence[x$correct_model],
          mean_coherence_est = x$mean_coherence[indicator],
          var_coherence_correct = x$var_coherence[x$correct_model],
          var_coherence_est = x$var_coherence[indicator],
          mean_prevalence_correct = x$mean_prevalence[x$correct_model],
          mean_prevalence_est = x$mean_prevalence[indicator],
          var_prevalence_correct = x$var_prevalence[x$correct_model],
          var_prevalence_est = x$var_prevalence[indicator],
          geweke_correct = x$geweke_stat[x$correct_model],
          geweke_est = x$geweke_stat[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()

# if eta is flat
eta_flat <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_flat_eta
      
      if (TRUE) {
        out <- tibble(
          id = x$id[1],
          ll_correct = x$ll[x$correct_model],
          ll_est = x$ll[indicator],
          r2_correct = x$r2[x$correct_model],
          r2_est = x$r2[indicator],
          mean_coherence_correct = x$mean_coherence[x$correct_model],
          mean_coherence_est = x$mean_coherence[indicator],
          var_coherence_correct = x$var_coherence[x$correct_model],
          var_coherence_est = x$var_coherence[indicator],
          mean_prevalence_correct = x$mean_prevalence[x$correct_model],
          mean_prevalence_est = x$mean_prevalence[indicator],
          var_prevalence_correct = x$var_prevalence[x$correct_model],
          var_prevalence_est = x$var_prevalence[indicator],
          geweke_correct = x$geweke_stat[x$correct_model],
          geweke_est = x$geweke_stat[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()

# if alpha_sum is too big
alpha_sum_big <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_alpha_sum != x$alpha_sum
      
      if (x$est_alpha_sum[x$correct_model] < x$est_alpha_sum[indicator]) {
        out <- tibble(
          id = x$id[1],
          ll_correct = x$ll[x$correct_model],
          ll_est = x$ll[indicator],
          r2_correct = x$r2[x$correct_model],
          r2_est = x$r2[indicator],
          mean_coherence_correct = x$mean_coherence[x$correct_model],
          mean_coherence_est = x$mean_coherence[indicator],
          var_coherence_correct = x$var_coherence[x$correct_model],
          var_coherence_est = x$var_coherence[indicator],
          mean_prevalence_correct = x$mean_prevalence[x$correct_model],
          mean_prevalence_est = x$mean_prevalence[indicator],
          var_prevalence_correct = x$var_prevalence[x$correct_model],
          var_prevalence_est = x$var_prevalence[indicator],
          geweke_correct = x$geweke_stat[x$correct_model],
          geweke_est = x$geweke_stat[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()

# if alpha_sum is too small
alpha_sum_small <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_alpha_sum != x$alpha_sum
      
      if (x$est_alpha_sum[x$correct_model] > x$est_alpha_sum[indicator]) {
        out <- tibble(
          id = x$id[1],
          ll_correct = x$ll[x$correct_model],
          ll_est = x$ll[indicator],
          r2_correct = x$r2[x$correct_model],
          r2_est = x$r2[indicator],
          mean_coherence_correct = x$mean_coherence[x$correct_model],
          mean_coherence_est = x$mean_coherence[indicator],
          var_coherence_correct = x$var_coherence[x$correct_model],
          var_coherence_est = x$var_coherence[indicator],
          mean_prevalence_correct = x$mean_prevalence[x$correct_model],
          mean_prevalence_est = x$mean_prevalence[indicator],
          var_prevalence_correct = x$var_prevalence[x$correct_model],
          var_prevalence_est = x$var_prevalence[indicator],
          geweke_correct = x$geweke_stat[x$correct_model],
          geweke_est = x$geweke_stat[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()

# calculate a bunch of t-statistics
t_helper <- function(correct, est) {
  
  my_t <- t.test(est, correct, paired = TRUE)
  
  tibble(
    mean_diff = my_t$estimate,
    ci_95_low = my_t$conf.int[1],
    ci_95_high = my_t$conf.int[2],
    t_stat =  my_t$statistic
  )
  
}

df_helper <- function(df, par) {
  
  out <- list(
    ll = t_helper(df$ll_correct, df$ll_est) |>
      mutate(outcome =  "ll"),
    r2 = t_helper(df$r2_correct, df$r2_est) |>
      mutate(outcome =  "r2"),
    mean_coherence = t_helper(df$mean_coherence_correct, df$mean_coherence_est) |>
      mutate(outcome =  "mean_coherence"),
    mean_prevalence = t_helper(df$mean_prevalence_correct, df$mean_prevalence_est) |>
      mutate(outcome =  "mean_prevalence"),
    geweke = t_helper(df$geweke_correct, df$geweke_est) |>
      mutate(outcome ="geweke")
  ) |>
    bind_rows() |>
    mutate(
      parameter = par
    )
  
  return(out)
}

t_tests <-
  list(
    alpha_sum_big = df_helper(alpha_sum_big, "alpha_sum_big"),
    alpha_sum_small = df_helper(alpha_sum_small, "alpha_sum_small"),
    eta_flat = df_helper(eta_flat, "eta_flat"),
    eta_sum_big = df_helper(eta_sum_big, "eta_sum_big"),
    eta_sum_small = df_helper(eta_sum_small, "eta_sum_small"),
    nk_big = df_helper(nk_big, "nk_big"),
    nk_small = df_helper(nk_small, "nk_small")
  ) |>
  bind_rows()

### write out object(s) for analysis ----

write_rds(
  x = t_tests,
  file = "data-derived/zipf-analysis/model-t-tests.rds"
)

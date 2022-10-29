# This script prepares/analyzes data comparing correct/incorrect model specifications


### load libraries ---
library(tidyverse)
library(tidylda)
library(Matrix)

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

### Augment metrics table with partial r2 and lambda coherence ----

# define some functions 
remove_one_topic <- function(model, topic_index) {
  
  new_theta <- model$theta[, - topic_index]
  
  new_theta <- new_theta / rowSums(new_theta, na.rm = TRUE)
  
  new_theta[is.na(new_theta)] <- 0
  
  new_beta <- model$beta[- topic_index, ]
  
  list(
    theta = new_theta,
    beta = new_beta
  )
  
}

get_partial_r2 <-
  function(m, dtm) {
    partial_r2 <- 
      tibble(
        topic = 1:nrow(m$beta),
        partial_r2 = (1:nrow(m$beta)) |> 
          map(function(x) {
            m_new <- remove_one_topic(m, x)
            
            r2 <- tidylda:::calc_lda_r2(
              dtm = dtm, 
              theta = m_new$theta,
              beta = m_new$beta,
              threads = parallel::detectCores() - 1
            )
          })
      )
    
    partial_r2 <- 
      partial_r2 |>
      mutate(
        partial_r2 = m$r2 - unlist(partial_r2)
      )
    
    partial_r2
  }

# get objects to loop
model_names <-
  list.files("data-derived/zipf-analysis/models")

augmentation <-
  model_names |>
  parallel::mclapply(
    function(filename) {
      
      # load model
      models <- read_rds(
        paste0("data-derived/zipf-analysis/models/", filename)
      )
      
      # load dtm
      dat <- read_rds(
        paste0("data-derived/zipf-analysis/simulated-data/", filename)
      )
      
      # calculate relevant stats
      out <-
        models |>
        map(
          function(m){
            
            get_partial_r2(m = m, dtm = dat$dtm[, colnames(m$beta)]) |>
              mutate(
                lambda_coherence = tidylda:::calc_prob_coherence(
                  beta = m$lambda,
                  data = dat$dtm[, colnames(m$beta)]
                ),
                coherence = m$summary$coherence,
                prevalence = m$summary$prevalence
              )
          }
        )
      
      # exit
      out
    }, mc.cores = parallel::detectCores() - 1
  )

write_rds(
  augmentation,
  file = "data-derived/zipf-analysis/partial-r2.rds"
)

model_metrics <- 
  model_metrics |>
  cbind(
    augmentation |>
      map(function(x){
        x |> map(function(y){
          y |> summarize(
            mean_partial_r2 = mean(partial_r2, na.rm = T),
            mean_lambda_coherence = mean(lambda_coherence, na.rm = T)
          )
        }) |> bind_rows()
      }) |> bind_rows()
  ) |> 
  as_tibble()




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
          geweke_est = x$geweke_stat[indicator],
          mean_partial_r2_correct = x$mean_partial_r2[x$correct_model],
          mean_partial_r2_est = x$mean_partial_r2[indicator],
          mean_lambda_coherence_correct = x$mean_lambda_coherence[x$correct_model],
          mean_lambda_coherence_est = x$mean_lambda_coherence[indicator]
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
          geweke_est = x$geweke_stat[indicator],
          mean_partial_r2_correct = x$mean_partial_r2[x$correct_model],
          mean_partial_r2_est = x$mean_partial_r2[indicator],
          mean_lambda_coherence_correct = x$mean_lambda_coherence[x$correct_model],
          mean_lambda_coherence_est = x$mean_lambda_coherence[indicator]
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
          geweke_est = x$geweke_stat[indicator],
          mean_partial_r2_correct = x$mean_partial_r2[x$correct_model],
          mean_partial_r2_est = x$mean_partial_r2[indicator],
          mean_lambda_coherence_correct = x$mean_lambda_coherence[x$correct_model],
          mean_lambda_coherence_est = x$mean_lambda_coherence[indicator]
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
          geweke_est = x$geweke_stat[indicator],
          mean_partial_r2_correct = x$mean_partial_r2[x$correct_model],
          mean_partial_r2_est = x$mean_partial_r2[indicator],
          mean_lambda_coherence_correct = x$mean_lambda_coherence[x$correct_model],
          mean_lambda_coherence_est = x$mean_lambda_coherence[indicator]
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
          geweke_est = x$geweke_stat[indicator],
          mean_partial_r2_correct = x$mean_partial_r2[x$correct_model],
          mean_partial_r2_est = x$mean_partial_r2[indicator],
          mean_lambda_coherence_correct = x$mean_lambda_coherence[x$correct_model],
          mean_lambda_coherence_est = x$mean_lambda_coherence[indicator]
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
          geweke_est = x$geweke_stat[indicator],
          mean_partial_r2_correct = x$mean_partial_r2[x$correct_model],
          mean_partial_r2_est = x$mean_partial_r2[indicator],
          mean_lambda_coherence_correct = x$mean_lambda_coherence[x$correct_model],
          mean_lambda_coherence_est = x$mean_lambda_coherence[indicator]
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
          geweke_est = x$geweke_stat[indicator],
          mean_partial_r2_correct = x$mean_partial_r2[x$correct_model],
          mean_partial_r2_est = x$mean_partial_r2[indicator],
          mean_lambda_coherence_correct = x$mean_lambda_coherence[x$correct_model],
          mean_lambda_coherence_est = x$mean_lambda_coherence[indicator]
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
      mutate(outcome ="geweke"),
    mean_partial_r2 = t_helper(df$mean_partial_r2_correct, df$mean_partial_r2_est) |>
      mutate(outcome ="partial_r2"),
    mean_lambda_coherence = t_helper(df$mean_lambda_coherence_correct, df$mean_lambda_coherence_est) |>
      mutate(outcome ="mean_lambda_coherence")
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

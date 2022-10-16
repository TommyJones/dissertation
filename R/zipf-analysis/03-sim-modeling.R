# this script creates several topic models for simulated data sets
# the purpose to is to see the effect of mis-specifying each hyperparameter
# while holding the others constant.
# Looking at the effects of:
#   Over/under estimating Nk
#   Over/under estimating sum_eta
#   Over/under estimating sum_alpha
#   Making sum_eta flat when it's a PL


set.seed(3851466)

### load libraries ----
library(tidyverse)
library(Matrix)
library(tidylda)
library(tmsamples)
library(coda)

### list of hyper parameters to sample from for misspecifications ----
sim_metrics <- 
  list(
    Nk = c(25, 50, 100, 200),
    eta_sum = c(50, 250, 500, 1000),
    alpha_sum = c(0.5, 1, 3, 5)
  )

### load list of dtms to loop over ----
dtm_list <- 
  list.files(
    path = "data-derived/zipf-analysis/simulated-data",
    full.names = TRUE
  )


### sample dtms so we get even distribution among parallel threads ----
dtm_list <- dtm_list[sample(seq_along(dtm_list), length(dtm_list))]


### divide DTMs into batches so we have fewer parallel jobs ----
# this lowers switching costs and reduces risks of side effects that have
# killed this job in the past

batch_size <-
  (length(dtm_list) / parallel::detectCores() - 1) |>
  floor()

batches <- seq(1, length(dtm_list), by = batch_size)

dtm_list <-
  batches |>
  map(
    function(x) {
      dtm_list[x:(min(x + batch_size - 1, length(dtm_list)))]
    }
  )


### loop over list of dtms and build models for each, saving metrics ----

metrics <- 
  dtm_list |>
  parallel::mclapply(
    function(batch) {
      
      out <-
        batch |>
        map(
          function(f) {
            d <- read_rds(f)
            
            # only keep dtm columns where a word was sampled at least once
            d$dtm <- d$dtm[, Matrix::colSums(d$dtm) > 0]
            
            
            # construct list of hyperparameters for each run
            # first is perfectly-specified
            # all others have some misspecification
            models <- 
              list(
                d$seed_pars |>
                  as_tibble() |>
                  mutate(flat_eta = FALSE) |>
                  select(
                    Nk,
                    eta_sum = beta_sum,
                    alpha_sum,
                    flat_eta
                  ),
                d$seed_pars |>
                  as_tibble() |>
                  mutate(flat_eta = FALSE) |>
                  mutate(
                    Nk = sample(sim_metrics$Nk[sim_metrics$Nk != Nk], 1)
                  ) |>
                  select(
                    Nk,
                    eta_sum = beta_sum,
                    alpha_sum,
                    flat_eta
                  ),
                d$seed_pars |>
                  as_tibble() |>
                  mutate(flat_eta = FALSE) |>
                  mutate(
                    beta_sum = sample(sim_metrics$eta_sum[sim_metrics$eta_sum != beta_sum], 1)
                  ) |>
                  select(
                    Nk,
                    eta_sum = beta_sum,
                    alpha_sum,
                    flat_eta
                  ),
                d$seed_pars |>
                  as_tibble() |>
                  mutate(flat_eta = FALSE) |>
                  mutate(
                    alpha_sum = sample(sim_metrics$alpha_sum[sim_metrics$alpha_sum != alpha_sum], 1)
                  ) |>
                  select(
                    Nk,
                    eta_sum = beta_sum,
                    alpha_sum,
                    flat_eta
                  ),
                d$seed_pars |>
                  as_tibble() |>
                  mutate(flat_eta = TRUE) |>
                  select(
                    Nk,
                    eta_sum = beta_sum,
                    alpha_sum,
                    flat_eta
                  )
              )
            
            
            # build models for each specification
            models <-
              models |>
              map(
                function(m) {
                  
                  # construct eta and alpha
                  if (! m$flat_eta) {
                    
                    eta <- generate_zipf(
                      vocab_size = ncol(d$dtm),
                      magnitude = m$eta_sum
                    )
                    
                  } else {
                    
                    eta <- m$eta_sum / ncol(d$dtm)
                    
                  }
                  
                  alpha <- m$alpha_sum / m$Nk
                  
                  
                  # build a model
                  lda <- 
                    tidylda(
                      data = d$dtm,
                      k = m$Nk,
                      iterations = 200, 
                      burnin = 150,
                      alpha = alpha,
                      eta = eta,
                      calc_likelihood = TRUE,
                      calc_r2 = TRUE
                    )
                  
                  # get geweke convergence diagnostic
                  g <- mcmc(lda$log_likelihood$log_likelihood[151:200]) |>
                    geweke.diag()
                  
                  # construct tibble of metrics
                  true_values <-
                    d$seed_pars |>
                    as_tibble()
                  
                  
                  model_values <- 
                    tibble(
                      est_k = m$Nk,
                      est_alpha_sum = m$alpha_sum,
                      est_eta_sum = m$eta_sum,
                      est_flat_eta = m$flat_eta,
                      r2 = lda$r2,
                      ll = mean(lda$log_likelihood$log_likelihood[151:200]),
                      geweke_stat = g$z,
                      mean_coherence = mean(lda$summary$coherence),
                      median_coherence = median(lda$summary$coherence),
                      var_coherence = var(lda$summary$coherence),
                      quantile_coherence_25 = quantile(lda$summary$coherence, .25),
                      quantile_coherence_75 = quantile(lda$summary$coherence, .75),
                      mean_prevalence = mean(lda$summary$prevalence),
                      median_prevalence = median(lda$summary$prevalence),
                      var_prevalence = var(lda$summary$prevalence),
                      quantile_prevalence_25 = quantile(lda$summary$prevalence, .25),
                      quantile_prevalence_75 = quantile(lda$summary$prevalence, .75)
                    )
                    
                  
                  metrics <- cbind(true_values, model_values)
                  
                  lda$metrics <- metrics
                  
                  lda
                }
              ) 
            
            # save the models
            write_rds(
              x = models,
              file = paste0(
                "data-derived/zipf-analysis/models/",
                paste(d$seed_pars, collapse = "-"),
                ".rds"
              ),
              compress = "gz"
            )
            
            # pull out the metrics, rbind them, then spit them out
            models |>
              map(function(m){
                m$metrics
              }) |>
              bind_rows()

          }
        )
      
      out
      
    },
    mc.cores = parallel::detectCores() - 1
  ) |> 
  bind_rows()

write_rds(
  x = metrics,
  file = "data-derived/zipf-analysis/model-metrics.rds"
)


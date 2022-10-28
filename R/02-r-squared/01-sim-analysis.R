## As of 2022-10-28, this script produces identical outputs to 
# R/zipf-analysis/01-sim-data.R So, I am archiving it.


# this script creates simulated data sets for the r-squared analysis
# convergence to a final model

library(tidyverse)
library(textmineR)
library(tmsamples)
library(furrr)

set.seed(8675201)

# set up futures for parallelism
# Edit: 2022-05-15 - using multithreading in C++ for parallelism b/c much faster
# plan(multisession, workers = future::availableCores() - 1, gc = TRUE)

### function to generate LDA parameters in a consistent way
gen_pars <- function(Nd, Nv, Nk, beta_sum, alpha_sum, alpha_shape = "flat") {
  
  z <- generate_zipf(vocab_size = Nv, magnitude = beta_sum)
  
  if (alpha_shape == "flat") {
    a <- rep(alpha_sum / Nk, Nk)
  } else {
    a <- rgamma(Nk, 1)
    
    a <- a / sum(a) * alpha_sum
  }
  
  par <- sample_parameters(alpha = a, beta = z, num_documents = Nd)
  
  par 
}

### Declare a function to calculate evaluation metrics ----
calc_metrics <- function(phi, theta, dtm, ...){
  K <- nrow(phi)
  D <- nrow(dtm)
  V <- ncol(dtm)
  len <- Matrix::rowSums(dtm)
  len <- mean(len)
  
  theta <- theta / rowSums(theta)
  phi <- phi / rowSums(phi)
  
  # these are for the "no model" ll in McFadden's
  phi_raw <- Matrix::colSums(dtm) + 1 / ncol(dtm)
  phi_raw <- phi_raw / sum(phi_raw)
  phi_raw <- rbind(phi_raw, phi_raw) 
  rownames(phi_raw) <- c("t_1", "t_2")
  colnames(phi_raw) <- colnames(dtm)
  
  theta_raw <- rep(0.5, nrow(dtm))
  theta_raw <- cbind(theta_raw, theta_raw)
  rownames(theta_raw) <- rownames(dtm)
  colnames(theta_raw) <- rownames(phi_raw)
  
  ll_model <- CalcLikelihood(dtm = dtm, phi = phi, theta = theta, ...)
  ll_raw <- CalcLikelihood(dtm = dtm, phi = phi_raw, theta = theta_raw, ...)
  r2 <- tidylda:::calc_lda_r2(dtm = dtm, theta = theta, beta = phi, threads = future::availableCores() - 1)
  r2_mac <- 1 - ll_model/ll_raw
  
  result <- list(
    ll_model = ll_model, 
    ll_raw = ll_raw, 
    r2 = r2, 
    r2_mac = r2_mac
  )
  
  return(result)
}

### Generate parameters, sample DTMs, get evaluation stats, and save to file ----
sim_metrics <- 
  list(
    Nv = c(1000, 5000, 10000, 20000),
    Nd = c(500, 1000, 2000, 4000),
    Nk = c(25, 50, 100, 200),
    doc_length = c(50, 100, 200, 400),
    beta_sum = c(50, 250, 500, 1000),
    alpha_sum = c(0.5, 1, 3, 5)
  ) |>
  cross() |>
  map(
    function(x) {
      
      # sample LDA parameters
      out <- list(
        seed_pars = x,
        lda_pars = gen_pars(
          Nd = x$Nd,
          Nv = x$Nv,
          Nk = x$Nk,
          beta_sum = x$beta_sum,
          alpha_sum = x$alpha_sum,
          alpha_shape = "flat"
        )
      )
      
      # sample DTM
      doc_lengths <- rpois(
        n = nrow(out$lda_pars$theta),
        lambda = out$seed_pars$doc_length
      )
      
      out$dtm <- try({
        sample_documents(
          theta = out$lda_pars$theta,
          phi = out$lda_pars$phi,
          doc_lengths = doc_lengths,
          verbose = FALSE,
          threads = future::availableCores() - 1
        )
      })
      
      # write output to file
      write_rds(
        out,
        file = paste0(
          "data-derived/r-squared/simulated-data/", 
          paste(out$seed_pars, collapse = "-"), 
          ".rds"
        )
      )
      
      # calculate evaluation metrics
      if (! inherits(out$dtm, "try-error")) {
        metrics <- 
          out$seed_pars |>
          as_tibble() |>
          cbind(
            calc_metrics(
              phi = out$lda_pars$phi,
              theta = out$lda_pars$theta,
              dtm = out$dtm
            ) |> as_tibble()
          )
      } else {
        metrics <-
          out$seed_pars |>
          as_tibble() |>
          cbind(
            tibble(
              ll_model = NA, 
              ll_raw = NA, 
              r2 = NA, 
              r2_mac = NA 
            )
          )
      }
      
      # spit out the metrics
      metrics
    }# , .options = furrr_options(seed = TRUE)
  ) |>
  bind_rows()

write_rds(
  sim_metrics,
  "data-derived/r-squared/sim-metrics.rds"
)




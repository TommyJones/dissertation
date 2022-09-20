# This script takes the DTMs of simulated data and calculates (potentially)
# relevant statistics. The goals of the data generated here are twofold:
#   1. To see which data generating specifications are consistent with Zipf's
#      and Heaps's laws (and possibly other properties of language)
#   2. To see if certain DTM statistics can predict / are associated with
#      certain hyperparameter choices

### load libraries ----
library(tidyverse)
library(Matrix)
library(poweRlaw)


### declare functions ----

# function uses regression to estimate heaps's law
# two problems: regression isn't a great way to do it
# and need to update how we get n and v based on sampling
estimate_heaps <- function(dtm) {
  
  heaps <- tibble(
    n = rowSums(dtm),
    v = rowSums(dtm > 0)
  )
  
  heaps <- lm(
    I(log(v)) ~ I(log(n)),
    data = heaps
  ) %>%
    summary() %>%
    .[["coefficients"]] %>%
    .[, "Estimate"]
  
  heaps
}

estimate_zipf <- function(dtm) {
  # estimate zipf coefficient from probability model
  
  word_freqs <- Matrix::colSums(dtm)
  
  m_m <- displ$new(word_freqs)
  
  est <- estimate_xmin(m_m)
  
  est
}

estimate_zipf2 <- function(dtm, xmin = 0) {
  # estimate zipf from a regression model
  
  word_freqs <- Matrix::colSums(dtm)
  
  word_freqs <- word_freqs[word_freqs >= xmin]
  
  est <- lm(
    I(log10(sort(word_freqs, decreasing = TRUE))) ~ I(log10(seq_along(word_freqs)))
  ) %>%
    summary() %>%
    .[["coefficients"]] %>%
    .[, "Estimate"]
  
  est
  
}


calc_dtm_metrics <- function(dtm) {
  
  vocab <- sum(Matrix::colSums(dtm) > 0)
  
  word_vol <- sum(dtm)
  
  num_docs <- nrow(dtm)
  
  word_var_in_docs <- numeric(num_docs)
  
  for (j in seq_len(num_docs)) {
    word_var_in_docs[j] <- var(dtm[j, ])
  }
  
  doc_lengths <- Matrix::rowSums(dtm)
  
  heaps <- 
    dtm |> 
    estimate_heaps()
  
  zipf <- 
    dtm |> 
    estimate_zipf()
  
  zipf2 <- 
    dtm |> 
    estimate_zipf2(xmin = zipf$xmin)
  
  
  dtm_stats <- 
    tibble(
      vocab = vocab,
      word_vol = word_vol,
      mean_doc_length = mean(doc_lengths),
      var_doc_length = var(doc_lengths),
      mean_word_var_in_docs = mean(word_var_in_docs),
      var_word_var_in_docs = var(word_var_in_docs),
      heaps_const = heaps[1],
      heaps_coef = heaps[2],
      zipf_xmin = zipf$xmin,
      zipf_coef = zipf$pars,
      zipf2_const = zipf2[1],
      zipf2_coef = zipf2[2]
    )
  
  dtm_stats
  
}

### get list of DTMs to load ----
file_list <- 
  list.files(
    "data-derived/zipf-analysis/simulated-data",
    full.names = TRUE
    )


### loop over the DTMs and calculate metrics ----

metrics <- 
  file_list |>
  parallel::mclapply(
    function(f) {
      
      d <- read_rds(f)
      
      seed_pars <- 
        d$seed_pars |>
        as_tibble()
      
      dtm_metrics <- 
        d$dtm |>
        calc_dtm_metrics()
      
      out <- cbind(seed_pars, dtm_metrics)
      
      out
      
    }, mc.cores = parallel::detectCores() - 1
  ) |>
  bind_rows()

write_rds(
  x = metrics,
  file = "data-derived/zipf-analysis/dtm-metrics.rds"
)




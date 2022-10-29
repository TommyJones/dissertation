# this script creates an NIH topic model for r-squared analysis

library(tidyverse)
library(textmineR)
library(tidytext)
library(tidylda)
library(furrr)

set.seed(3851466)

# load nih data
nih <- read_csv("data-raw/RePORTER_PRJABS_C_FY2014.csv")

names(nih) <- tolower(names(nih))

# sample 1000 documents and create a dtm
nih_dtm <- 
  nih[sample(seq_len(nrow(nih)), 1000), ] |>
  unnest_tokens(
    output = word, 
    input = abstract_text,
    stopwords = stop_words$word,
    token = "words"
  ) |> 
  filter(! is.na(word)) |>
  count(application_id, word) |>
  filter(n>1) |> #Filtering for words/bigrams per document, rather than per corpus
  cast_sparse(application_id, word, n)


write_rds(
  nih_dtm,
  "data-derived/r-squared/nih-dtm.rds"
)


# build models from k = 50 to k = 500
k_list <- seq(from = 50, to=500, by=50)

# plan(
#   multisession, 
#   workers = min(length(k_list), future::availableCores() - 1), 
#   gc = TRUE
# ) # set up futures for parallelism


nih_models <- 
  k_list |>
  parallel::mclapply(
    function(k) {
      
      m <- tidylda(
        data = nih_dtm,
        k = k,
        iterations = 200,
        burnin = 150,
        calc_r2 = TRUE
      )
      
    }, mc.cores = min(length(k_list), future::availableCores() - 1) # options = furrr_options(seed = TRUE)
  )

write_rds(
  nih_models,
  "data-derived/r-squared/nih-models.rds"
)

# pull out relevant stats into tibble
nih_metrics <- 
  nih_models |>
  map(
    function(x) {
      out <- 
        tibble(
          Nk = nrow(x$beta),
          ll = mean(x$log_likelihood$log_likelihood[151:200], na.rm = TRUE),
          r2 = x$r2
        )
      
      out
    }
  ) |> 
  bind_rows()

# add "naive" log likelihood
nih_metrics <- 
  nih_metrics |> 
  mutate(
    ll2 = nih_models |> 
      map(function(x) CalcLikelihood(nih_dtm, x$beta, x$theta)) |>
      unlist()
  )

# save nih-metrics
write_rds(
  nih_metrics,
  "data-derived/r-squared/nih-metrics.rds"
)


### partial R^2 analysis ----

# choose model with 100 topics

m <- nih_models[sapply(nih_models, function(x) nrow(x$beta) == 100)][[1]]

# define a function to remove one topic and rebalance
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

# loop over topics and get r-squared
partial_r2 <- 
  tibble(
    topic = 1:nrow(m$beta),
    partial_r2 = (1:nrow(m$beta)) |> 
      map(function(x) {
        m_new <- remove_one_topic(m, x)
        
        r2 <- tidylda:::calc_lda_r2(
          dtm = nih_dtm, 
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

write_rds(
  partial_r2,
  "data-derived/r-squared/nih-partial-r2.rds"
)


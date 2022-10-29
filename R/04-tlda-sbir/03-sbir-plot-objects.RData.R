# This script constructs objects for plotting from the sbir time series models

# load libraries
library(tidyverse)
library(textmineR)
library(tidylda)
library(tidytext)

# load("data-derived/tlda-sbir/sbir_models-by-award-year.RData")

sbir_models <- read_rds("data-derived/tlda-sbir/models-a-0.8.rds")

# get prevalence by year from the summary tables ----
year_prev <- 
  sbir_models %>%
  map(function(x){
    out <- x$summary
  })

for (j in seq_along(year_prev)) {
  year_prev[[j]] <- 
    year_prev[[j]] %>%
    mutate(
      year = as.numeric(names(year_prev)[j])
    )
}

year_prev <- 
  year_prev %>%
  bind_rows()


# get top words for each topic over time ----
year_words <- 
  sbir_models %>%
  parallel::mclapply(function(x) {
    
    # Beta
    beta <- tidy(x, matrix = "beta")
    
    # old school R here
    beta <- by(
      beta, 
      INDICES = beta$topic, 
      function(y) {
        y %>%
          arrange(desc(beta)) %>%
          .[1:10, ]
      })
    
    beta <- do.call(rbind, beta)
    
    # Lambda
    lambda <- tidy(x, matrix = "lambda")
    
    # old school R here
    lambda <- by(
      lambda, 
      INDICES = lambda$topic, 
      function(y) {
        y %>%
          arrange(desc(lambda)) %>%
          .[1:10, ]
      })
    
    lambda <- do.call(rbind, lambda)
    
    out <- beta %>%
      full_join(
        lambda
      )
    
  }, mc.cores = 7)

for (j in seq_along(year_words)) {
  year_words[[j]] <- 
    year_words[[j]] %>%
    mutate(year = as.numeric(names(year_words)[j]))
}

year_words <- 
  year_words %>%
  bind_rows()


# Get aggregate evaluation statistics over time ----

agg_stats <- tibble(
  year = as.numeric(names(sbir_models)),
  r2 = sbir_models |> sapply(function(x) x$r2),
  coherence_weighted = sbir_models |> 
    sapply(function(x) {
      sum(x$summary$coherence * x$summary$prevalence, na.rm = TRUE) / 
        nrow(x$summary)
    }) ,
  coherence_unweighted = sbir_models |> 
    sapply(function(x) mean(x$summary$coherence, na.rm = TRUE))
)

# Get hellinger distance changes, year to year and topic to topic ----

# define a function to get the change in hellinger distance from one
# model to the next

get_hellinger_delta <- function(model1, model2) {
  
  vocab <- intersect(
    colnames(model1$beta),
    colnames(model2$beta)
  )
  
  out <- sapply(1:nrow(model1$beta), function(x) {
    textmineR::CalcHellingerDist(
      x = model1$beta[x, vocab],
      y = model2$beta[x, vocab]
    )
  }) 
  
  out
  
}

# calculate hellinger distance for all
hellinger <- vector(mode = "list", length = length(sbir_models))

hellinger[[1]] <- tibble(
  year = rep(as.numeric(names(sbir_models)[1]), nrow(sbir_models[[1]]$beta)),
  topic = as.numeric(rownames(sbir_models[[1]]$beta)),
  hellinger = NA
)

for (j in 2:length(hellinger)) {
  
  hellinger[[j]] <- tibble(
    year = rep(as.numeric(names(sbir_models)[j]), nrow(sbir_models[[j - 1]]$beta)),
    topic = as.numeric(rownames(sbir_models[[j - 1]]$beta)),
    hellinger = get_hellinger_delta(
      model1 = sbir_models[[j - 1]],
      model2 = sbir_models[[j]]
    )
  )
  
}

hellinger <- 
  hellinger |>
  bind_rows()

# save outcomes ----
save(
  year_prev,
  year_words,
  agg_stats,
  hellinger,
  file = "data-derived/tlda-sbir/plot-objs.RData"
)

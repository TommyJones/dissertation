# This script builds topic models over time for the SBIR corpus
# It varries a over {0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6}
# smaller range than on the simulations because we know that 0.8 is kind of the
# critical area based on the simulation analyses

set.seed(42)

# load libraries
library(tidyverse)
library(textmineR)
library(tidylda)
library(tidytext)

# load cleaned data
load("data-derived/tlda-sbir/sbir-clean.RData")

# divide by time
# build a model annually
# we can still get monthly (or daily) assignments, the models just don't update
# that frequently

sbir_text <- by(
  data = sbir_text,
  INDICES = sbir$award_year,
  FUN = function(x) {
    x
  },
  simplify = FALSE
)

# DTM it all
sbir_dtms <- 
  sbir_text %>%
  parallel::mclapply(function(x){
    
    # tokenize using tidytext's unnest_tokens
    tidy_docs <- x %>% 
      mutate(text = paste(award_title, abstract)) %>%
      select(sbir_id, text) %>%
      unnest_tokens(output = word, 
                    input = text,
                    stopwords = stop_words$word,
                    token = "words") %>% 
      filter(! is.na(word)) %>%
      count(sbir_id, word) %>% 
      filter(n>1) #Filtering for words/bigrams per document, rather than per corpus
    
    dtm <- 
      tidy_docs %>%
      cast_sparse(sbir_id, word, n)
    
    dtm
    
  },
  mc.cores = 4)

save(
  sbir_dtms, 
  file = "data-derived/tlda-sbir/sbir-dtm-by-year.RData"
)


# fit initial 1983 model

m_1983 <- tidylda(
  data = sbir_dtms[[1]],
  k = 100,
  calc_likelihood = TRUE,
  calc_r2 = TRUE,
  iterations = 200,
  burnin = 150,
  verbose = TRUE
)


# definie a function to get the change in hellinger distance from one
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


# initialize a list to hold models
a_range <- seq(0.4, 1.6, by = 0.2)

names(a_range) <- a_range


# initialize a list of iterators to run in parallel
model_metrics <- map(
  a_range,
  function(a) {
    
    list(
      a = a,
      m_1983 = m_1983,
      sbir_dtms = sbir_dtms
    )
  }
)

# loop to build multiple model chains
# should be parallel, but for some unknown reason processes die using mclapply
model_metrics <- map( 
  model_metrics,
  function(m) {
    
    sbir_models <- vector(mode = "list", length = length(m$sbir_dtms))
    
    names(sbir_models) <- names(m$sbir_dtms)
    
    sbir_models[[1]] <- m$m_1983
    
    
    for (j in 2:length(sbir_models)) {
      sbir_models[[j]] <- refit(
        object = sbir_models[[j - 1]],
        new_data = m$sbir_dtms[[j]],
        iterations = 200,
        burnin = 150,
        prior_weight = m$a,
        calc_likelihood = TRUE,
        calc_r2 = TRUE,
        additional_k = 1,
        verbose = TRUE
      )
      cat("\n", names(sbir_models)[j], "\n")
    }
    
    write_rds(
      x = sbir_models,
      file = paste0("data-derived/tlda-sbir/models-a-",m$a, ".rds")
    )
    
    # gc()
    
    # get metrics from that model
    model_metrics <-
      map(
        .x = sbir_models,
        .f = function(.x) {
          
          out <- tibble(
            a = m$a,
            r2 = .x$r2,
            likelihood = mean(.x$log_likelihood$log_likelihood[51:200], na.rm = TRUE),
            mean_coh = mean(.x$summary$coherence, na.rm = TRUE),
            var_coh = var(.x$summary$coherence, na.rm = TRUE)
          )
          
          out
          
        }
      ) |>
      bind_rows()
    
    model_metrics <- 
      model_metrics |>
      mutate(
        year = sbir_models |> names() |> as.numeric(),
        mean_hellinger = NA,
        var_hellinger = NA
      )
    
    for (j in 2:length(sbir_models)) {
      
      h <- get_hellinger_delta(
        model1 = sbir_models[[j - 1]],
        model2 = sbir_models[[j]]
      )
      
      model_metrics$mean_hellinger[j] <- mean(h, na.rm = TRUE)
      
      model_metrics$var_hellinger[j] <- var(h, na.rm = TRUE)
      
    }
    
    model_metrics
    
  }
)


model_metrics <-
  model_metrics |> 
  bind_rows()

write_rds(
  model_metrics,
  "data-derived/tlda-sbir/model-metrics.rds"
)

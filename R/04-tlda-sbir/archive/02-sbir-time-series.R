# This script builds topic models over time for the SBIR corpus

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

# fit topic models for each year, with the last year's model as prior
sbir_models <- vector(mode = "list", length = length(sbir_dtms))

names(sbir_models) <- names(sbir_dtms) # name is year ending

sbir_models[[1]] <- tidylda(
  data = sbir_dtms[[1]],
  k = 100,
  calc_likelihood = TRUE,
  calc_r2 = TRUE,
  iterations = 200,
  burnin = 150,
  verbose = TRUE
)

for (j in 2:length(sbir_models)) {
  sbir_models[[j]] <- refit(
    object = sbir_models[[j - 1]],
    new_data = sbir_dtms[[j]],
    iterations = 200,
    burnin = 150,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    additional_k = 2,
    verbose = TRUE
  )
  cat("\n", j, "\n")
}

save(
  sbir_models,
  file = "data-derived/tlda-sbir/sbir_models-by-award-year.RData"
)

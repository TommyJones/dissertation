# this script is some experimental code for calculating descriptive stats from
# DTMs/TCMs

library(tidyverse)
library(Matrix)
library(poweRlaw)

# using NIH sample DTM for this experiment
dtm <- textmineR::nih_sample_dtm

# create some stats

vocab <- ncol(dtm)

word_vol <- sum(dtm)

num_docs <- nrow(dtm)

word_var_in_docs <- numeric(num_docs)

for (j in seq_len(num_docs)) {
  word_var_in_docs[j] <- var(dtm[j, ])
}

doc_lengths <- Matrix::rowSums(dtm)


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

heaps <- dtm %>% estimate_heaps()

estimate_zipf <- function(dtm) {
  # estimate zipf coefficient from probability model
  
  word_freqs <- Matrix::colSums(dtm)
  
  m_m <- displ$new(word_freqs)
  
  est <- estimate_xmin(m_m)
  
  est
}

zipf <- dtm %>% estimate_zipf()

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

zipf2 <- dtm %>% estimate_zipf2(xmin = zipf$xmin)

# create an object to store descriptive stats
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



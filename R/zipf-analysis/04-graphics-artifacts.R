library(tidyverse)
library(Matrix)
library(tmsamples)
library(tidytext)

set.seed(8675309)

### create a function for simulating a zipf distribution from a corpus ----
# function to imitate zipf distribution of a corpus
imitate_zipf <- function(dtm, sum_eta = 1000) {
  
  # get word freqs
  wf <- Matrix::colSums(dtm)
  
  # get xmin
  m <- displ$new(wf)
  
  xmin <- estimate_xmin(m)
  
  # do log-log regression
  d <- tibble(
    rank = seq_along(wf),
    freq = sort(wf, decreasing = TRUE)
  ) 
  
  f <- lm(
    log(freq) ~ I(log(rank)),
    data = d |>
      filter(freq >= xmin$xmin)
  )
  
  # solve for eta_v / sum(eta)
  d$p <- predict(f, newdata = d |> select(rank))
  
  eta <- exp(d$p - log(sum(dtm)))
  
  # make sure it sums to one
  eta <- eta / sum(eta)
  
  # scale by sum_eta
  eta <- eta * sum_eta
  
  # return from function
  eta
}


### create a matrix for plotting NIH vs simulated data ----
nih <- read_csv("data-raw/RePORTER_PRJABS_C_FY2014.csv.zip")

names(nih) <- tolower(names(nih))

# sample 1000 documents and create a dtm
nih_dtm <- 
  nih[sample(seq_len(nrow(nih)), 1000), ] |>
  unnest_tokens(
    output = word, 
    input = abstract_text,
    token = "ngrams",
    n_min = 1,
    n = 2
  ) |> 
  filter(! is.na(word)) |>
  count(application_id, word) |>
  filter(n>1) |> #Filtering for words/bigrams per document, rather than per corpus
  cast_sparse(application_id, word, n)

# get word frequencies
wf <- colSums(nih_dtm)

# create word frequencies of simulated data of same dimensions
pars1 <- sample_parameters(
  alpha = rep(0.1, 25),
  beta = imitate_zipf(
    dtm = nih_dtm,
    sum_eta = 1000
  ),
  num_documents = nrow(nih_dtm)
)

pars2 <- sample_parameters(
  alpha = generate_zipf( # from tm samples
    vocab_size = 25,
    magnitude = 0.1 * 25,
    zipf_par = 1
  ),
  beta = rep(0.05, ncol(nih_dtm)),
  num_documents = nrow(nih_dtm)
)

pars3 <- sample_parameters(
  alpha = generate_zipf(
    vocab_size = 25,
    magnitude = 0.1 * 25,
    zipf_par = 1
  ),
  beta = imitate_zipf(
    dtm = nih_dtm,
    sum_eta = 1000
  ),
  num_documents = nrow(nih_dtm)
)

pars4 <- sample_parameters(
  alpha = rep(0.1, 25),
  beta = rep(0.05, ncol(nih_dtm)),
  num_documents = nrow(nih_dtm)
)

# put together in a tibble for plotting
nih_plotmat <- 
  tibble(
    rank = seq_len(ncol(nih_dtm)),
    nih = nih_dtm |>
      colSums() |>
      sort(decreasing = TRUE),
    npl_pl = sample_documents(
      theta = pars1$theta,
      phi = pars1$phi,
      doc_lengths = rowSums(nih_dtm),
      threads = parallel::detectCores() - 1
    ) |>
      colSums() |>
      sort(decreasing = TRUE),
    pl_npl = sample_documents(
      theta = pars2$theta,
      phi = pars2$phi,
      doc_lengths = rowSums(nih_dtm),
      threads = parallel::detectCores() - 1
    ) |>
      colSums() |>
      sort(decreasing = TRUE),
    pl_pl = sample_documents(
      theta = pars3$theta,
      phi = pars3$phi,
      doc_lengths = rowSums(nih_dtm),
      threads = parallel::detectCores() - 1
    ) |>
      colSums() |>
      sort(decreasing = TRUE),
    npl_npl = sample_documents(
      theta = pars4$theta,
      phi = pars4$phi,
      doc_lengths = rowSums(nih_dtm),
      threads = parallel::detectCores() - 1
    ) |>
      colSums() |>
      sort(decreasing = TRUE)
  )

# example plot for me to copy to the rmd
# nih_plotmat |> ggplot() + 
#   geom_line(aes(x = log10(rank), y = log10(nih)), color = "red", lwd = 1.25) + 
#   geom_line(aes(x = log10(rank), y = log10(pl)), color = "blue", lty = 2, lwd = 1.25) +
#   geom_line(aes(x = log10(rank), y = log10(npl)), color = "green", lty = 3, lwd = 1.25)


# save results
write_rds(
  x = nih_plotmat,
  file = "data-derived/zipf-analysis/nih-plotmat.rds"
)




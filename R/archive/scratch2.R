library(Matrix)
library(poweRlaw)
library(tidyverse)
library(tidytext)
library(tmsamples)

set.seed(8675309)


### create a matrix for plotting NIH vs simulated data ----
nih <- read_csv("data-raw/RePORTER_PRJABS_C_FY2014.csv.zip")

names(nih) <- tolower(names(nih))

# sample 1000 documents and create a dtm
dtm <- 
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



# create word frequencies of simulated data of same dimensions
pars1 <- sample_parameters(
  alpha = rep(0.1, 25),
  beta = imitate_zipf(dtm, 1000),
  num_documents = nrow(dtm)
)

# put together in a tibble for plotting
nih_plotmat <- 
  tibble(
    rank = seq_len(ncol(dtm)),
    nih = dtm |>
      colSums() |>
      sort(decreasing = TRUE),
    npl_pl = sample_documents(
      theta = pars1$theta,
      phi = pars1$phi,
      doc_lengths = rowSums(dtm),
      threads = parallel::detectCores() - 1
    ) |>
      colSums() |>
      sort(decreasing = TRUE)
  )


nih_plotmat |> ggplot() + 
  geom_line(aes(x = rank, y = nih)) + 
  geom_line(aes(x = rank, y = npl_pl), color = "red") +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') 


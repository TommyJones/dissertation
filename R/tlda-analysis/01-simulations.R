# this script fits models across a range of "a" to see how it affects
# convergence to a final model

library(tidyverse)
library(tidytext)
library(tidylda)
library(tmsamples)

set.seed(8675201)


### create a common date/time for saving any elements to disk ----
clockmark <-
  Sys.time() %>%
  str_replace_all("[^a-zA-Z0-9]+", "-")


### generate population parameters ----
Nv = 5000

Nd = 10000

Nk = 25

gen_pars <- function(Nd, Nv, Nk, beta_sum, alpha_sum, alpha_shape = "flat") {
  
  z <- generate_zipf(vocab_size = Nv, magnitude = beta_sum)
  
  if (alpha_shape == "flat") {
    a <- rep(alpha_sum / Nk, Nk)
  } else {
    a <- rgamma(Nk, 1)
    
    a <- a / sum(a) * alpha_sum
  }
  
  par <- sample_parameters(alpha = a, beta = z, num_documents = Nd)
  
  par <- list(
    par = par,
    beta_sum = beta_sum,
    alpha_sum = alpha_sum,
    alpha_shape = alpha_shape
  )
}



beta_sums <- c(50, 250, 500, 1000)

alpha_sums <- c(0.5, 1, 3, 5)

alpha_shapes <- c("flat", "flat") # c("flat", "asymmetric")

doc_lengths <- c(50, 100, 200, 400)

pars <- cross3(
  alpha_sums, 
  beta_sums, 
  alpha_shapes
) %>%
  map(function(x) {
    names(x) <- c("alpha_sum", "beta_sum", "alpha_shape")
    x
  }) %>%
  map(function(x){
    gen_pars(
      Nd = Nd,
      Nv = Nv,
      Nk = Nk,
      beta_sum = x$beta_sum,
      alpha_sum = x$alpha_sum,
      alpha_shape = x$alpha_shape
    )
  }) %>%
  cross2(doc_lengths) %>%
  map(function(x){
    names(x) <- c("par", "doc_length")
    
    x
  })

# save pars to track results later
new_dir <- 
  paste0("data-derived/", clockmark, "/")

if (! dir.exists(new_dir))
  dir.create(new_dir)

save(
  pars,
  file = paste0(new_dir, "pop-pars.RData")
)

### sample document term matrices ----
dtms <- pars %>%
  parallel::mclapply(
    function(p) {
      # pull out relevant objects
      pop_pars <- p$par$par
      
      # generate documents
      doc_lengths <- rpois(
        n = nrow(pop_pars$theta),
        lambda = p$doc_length
      )
      
      docs <- sample_documents(
        theta = pop_pars$theta,
        phi = pop_pars$phi,
        doc_lengths = doc_lengths,
        threads = 1
      )
      
      # divide docs into batches
      batch_size <- 100
      
      batches <- seq(1, nrow(docs), by = batch_size)
      
      batches <- map(batches, function(x) x:min(x + batch_size - 1, nrow(docs))) 
      
      # list of dtms
      out <- batches %>%
        map(function(x) docs[x, ])
      
      out
    },
    mc.cores = parallel::detectCores() - 1
  )

save(
  dtms,
  file = paste0(new_dir, "sim-dtms.RData")
)

### across a range of "a", sample data and  fit models ----
# at each step, randomly shuffle the set of sampled DTMs
# this inserts randomness to guard against results being overfit on one set

a_range <- seq(0.2, 2, by = 0.2)

names(a_range) <- paste0("a-", a_range)

# add a number to pars to keep track of models downstream

# add model number to track results later
model_numbers <- seq_along(pars) %>% 
  as.character() %>%
  map(function(j) {
    if (nchar(j) == 1) {
      paste0("00", j)
    } else if (nchar(j) == 2) {
      paste0("0", j)
    } else {
      as.character(j)
    }
  })

for (j in seq_along(pars)) {
  pars[[j]]$model_num <- model_numbers[[j]]
}

# add dtms to pars
for (j in seq_along(pars)) {
  pars[[j]]$dtm <- dtms[[j]]
}

# re-save pars
save(
  pars,
  file = paste0(new_dir, "pop-pars.RData")
)

# calc Hellinger for real now
hdist <-
  a_range %>%
  map(
    function(a) {

      # helper function
      prep_dtm <- function(dtm) {
        # removes any vocabulary words that either appear in almost every document
        # or appear in no documents
        tf <- textmineR::TermDocFreq(dtm)
        
        vocab <- tf %>%
          filter(doc_freq > 0 & doc_freq < nrow(dtm) - 5) %>%
          select(term) %>%
          .[[1]]
        
        dtm[, vocab]
      }
      
      # declare hyper parameters
      k <- Nk + 5 # more estimated topics than actual
      
      alpha <- 0.05 # symmetric and small alpha no matter what
      
      eta <- 0.01 # symmetric and small eta no matter what
      
      # begin the modeling loop
      hdist <- 
        parallel::mcmapply(
          function(d_list, par_list) {
            
            models <- vector(mode = "list", length = length(d_list))
            
            models[[1]] <- try({
              tidylda(
                data = prep_dtm(d_list[[1]]),
                k = k,
                iterations = 100,
                burnin = 75,
                alpha = alpha,
                eta = eta,
                calc_likelihood = TRUE,
                calc_r2 = TRUE,
                verbose = TRUE
              )
            })
            
            for (j in 2:length(models)) {
              
              models[[j]] <- try({
                refit(
                  object = models[[j - 1]],
                  new_data = prep_dtm(d_list[[j]]),
                  iterations = 200,
                  burnin = 150,
                  prior_weight = a,
                  calc_likelihood = TRUE
                )
              })
              
            }
            
            # save models
            if (! exists(paste0(new_dir, "a-", a))) {
              dir.create(paste0(new_dir, "a-", a))
            }
            
            save(
              models,
              file = paste0(new_dir, "a-", a, "/m-", par_list$model_num, ".RData")
              )
            
            # calculate hellinger
            pop_pars <- par_list$par$par
            
            hdist <- models %>%
              map(function(m){
                try({
                  vocab <- intersect(colnames(m$beta), colnames(pop_pars$phi))
                  
                  h <- apply(m$beta[, vocab], 1, function(p){
                    apply(pop_pars$phi[, vocab], 1, function(q){
                      textmineR::CalcHellingerDist(p, q)
                    })
                  })
                  h
                })
              })
            
            # return hdist object
            hdist
            
          },
          mc.cores = parallel::detectCores() - 1,
          SIMPLIFY = FALSE,
          d_list = dtms,
          par_list = pars
        )
      
      hdist
    }
  )

save(
  hdist,
  file = paste0(new_dir, "hdist.RData")
)


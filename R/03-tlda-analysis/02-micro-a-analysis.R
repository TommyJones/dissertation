### This script runs models for a in 0.6 - 1 range ----


# initial setup
library(tidyverse)
library(tidytext)
library(tidylda)
library(tmsamples)

set.seed(8675201)

Nk <- 25

# load sim dtms and population parameters
pars <- read_rds("data-derived/tlda-sims/pop-pars.rds")

dtms <- read_rds("data-derived/tlda-sims/sim-dtms.rds")

# calc Hellinger for a smaller range of a values
a_range <- seq(0.65, 0.95, by = 0.05)

a_range <- a_range[a_range != 0.8]

names(a_range) <- paste0("a-", a_range)


cat(format(Sys.time(), "%a %b %d %X %Y"), " calculating Hellinger \n")

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
            if (! exists(paste0("data-derived/tlda-sims/", "a-", a))) {
              dir.create(paste0("data-derived/tlda-sims/", "a-", a))
            }
            
            save(
              models,
              file = paste0("data-derived/tlda-sims/", "a-", a, "/m-", par_list$model_num, ".RData")
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

# read in old hdist, append new hdist, and then save

hdist1 <- read_rds("data-derived/tlda-sims/hdist.rds")

hdist <- c(hdist1, hdist)

hdist <- hdist[order(names(hdist))]

write_rds(
  hdist,
  file = "data-derived/tlda-sims/hdist.rds"
)
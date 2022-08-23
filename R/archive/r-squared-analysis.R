##########################################################
# This script is, as of this writing, a holdover from the 
# self-contained R2 paper. I need to update its contents
# to work with the dissertation as a whole.
##########################################################

library(textmineR)

set.seed(8675309)

##########################################################
# Functions to generate a corpus with varying properties
###########################################################

### Function to simulate phi and theta parameters ------------------------------
SimulatePar <- function(D, V, K, alpha = rep(5 / K, K), beta_sum = 0.01 * V){
  
  # make a beta that follows zipf's law
  Zipf <- function(V) 1/(1:V * log(1.78 * V) ) # approximate empirical zipf distribution http://mathworld.wolfram.com/ZipfsLaw.html
  
  zipf.law <- Zipf(V=V)
  
  zipf.law <- zipf.law / sum(zipf.law) # normalize to sum to one (rounding error makes this slightly larger than one)
  
  beta <- zipf.law * beta_sum
  
  
  
  # phi
  
  # Use dirichlet, not multinomial. Variability of words within documents is
  # controlled by scale of beta. If smallest entries of beta are tiny, then
  # you have a "small sample" and zipf's law will not hold
  phi <- gtools::rdirichlet(n = K, alpha = beta)
  
  
  rownames(phi) <- paste("t", 1:K, sep="_")
  
  if (is.null(names(beta))){
    colnames(phi) <- paste("w", 1:V, sep = "_")
  } else {
    colnames(phi) <- names(beta)
  }
  
  # theta
  
  theta <- gtools::rdirichlet(n = D, alpha = alpha)
  
  # theta <- t(theta)
  rownames(theta) <- paste("d", 1:D, sep="_")
  colnames(theta) <- rownames(phi)
  
  return(list(theta=theta, phi=phi))
}

### Function to sample from phi and theta to construct a dtm -------------------
SampleDocs <- function(phi, theta, lambda = NULL, cpus = 4){
  
  D <- nrow(theta)
  K <- ncol(theta)
  V <- ncol(phi)
  
  # vector of document lengths
  if ( ! is.null(lambda) & length(lambda) == 1) {
    
    doc_lengths <- rpois(n = D, lambda = lambda) # document length / number of samples
    
  } else if (! is.null(lambda) & length(lambda) > 1) {
    doc_lengths <- lambda
  } else {
    
    warning("ill-defined lambda specified, choosing lambda = 300")
    
    doc_lengths <- rpois(n = 1, lambda = 300) # document length / number of samples
    
  }
  
  # make theta into an iterator
  docnames <- rownames(theta)
  
  step <- ceiling(nrow(theta) / cpus)
  
  batches <- seq(1, nrow(theta), by = step)
  
  theta <- lapply(batches, function(x){
    theta[ x:min(x + step - 1, nrow(theta)) , ]
  })
  
  theta <- lapply(theta, function(x){
    as.list(as.data.frame(t(x)))
  })
  
  # iterate over theta to perform sampling of words
  dtm <- parallel::mclapply(theta, function(batch){
    
    result <- lapply(batch, function(d){
      
      n_d <- sample(doc_lengths, 1)
      
      # take n_d topic samples
      topics <- rmultinom(n = n_d, size = 1, prob = d)
      
      # reduce and format to optimize for time/memory
      topics <- rowSums(topics)
      
      topics <- rbind(topic_index = seq_along(topics), 
                      times_sampled = topics)
      
      topics <- topics[ , topics[ "times_sampled" , ] > 0 ]
      
      # sample words from each topic
      
      if (is.null(dim(topics))) # corrects for event only one topic is sampled
        x <- matrix(x, nrow = 2)
      
      words <- apply(topics, 2, function(x){
        
        result <- rmultinom(n = x[ "times_sampled" ], size = 1, prob = phi[ x[ "topic_index" ] , ])
        
        if (! is.null(dim(result)))
          result <- rowSums(result)
        
        result
      })
      
      # turn into a single vector
      words <- Matrix::Matrix(rowSums(words), nrow = 1, sparse = T)
      
      # garbage collection to reduce memory footprint
      gc()
      
      return(words)
    })
    
    result <- do.call(rbind, result)
  }, mc.cores = cpus)
  
  dtm <- do.call(rbind, dtm)
  
  rownames(dtm) <- docnames
  colnames(dtm) <- colnames(phi)
  
  return(dtm)
}

### Function to calculate metrics ----
MetricFun <- function(phi, theta, dtm, ...){
  K <- nrow(phi)
  D <- nrow(dtm)
  V <- ncol(dtm)
  len <- Matrix::rowSums(dtm)
  len <- mean(len)
  
  theta <- theta / rowSums(theta)
  phi <- phi / rowSums(phi)
  
  # these are for the "no model" ll in McFadden's
  phi_raw <- Matrix::colSums(dtm) + 1 / ncol(dtm)
  phi_raw <- phi_raw / sum(phi_raw)
  phi_raw <- rbind(phi_raw, phi_raw) 
  rownames(phi_raw) <- c("t_1", "t_2")
  colnames(phi_raw) <- colnames(dtm)
  
  theta_raw <- rep(0.5, nrow(dtm))
  theta_raw <- cbind(theta_raw, theta_raw)
  rownames(theta_raw) <- rownames(dtm)
  colnames(theta_raw) <- rownames(phi_raw)
  
  ll_model <- CalcLikelihood(dtm=dtm, phi=phi, theta=theta, ...)
  ll_raw <- CalcLikelihood(dtm=dtm, phi=phi_raw, theta=theta_raw, ...)
  r2 <- CalcTopicModelR2(dtm=dtm, phi=phi, theta=theta, ...)
  r2_mac <- 1 - ll_model/ll_raw
  
  result <- data.frame(K=K, D=D, V=V, len=len, ll_model=ll_model, ll_raw=ll_raw, r2=r2, r2_mac=r2_mac, stringsAsFactors=FALSE)
  return(result)
}

##########################################################
# Actual executable code and analysis
##########################################################

### Default parameters ----
K <- 50
D <- 2000
V <- 5000
lambda <- 500


### Vary K ----
# lists of parameters to vary
k_list <- c(20, 50, 100, 200)
names(k_list) <- paste("k", k_list, sep="_")


par <- parallel::mclapply(k_list, 
                          function(k) SimulatePar(D = D, V = V, K = k),
                          mc.cores = parallel::detectCores() - 1)

k_m <- lapply(par, function(p){
  
  d <- SampleDocs(phi = p$phi, theta = p$theta, lambda = lambda)
  
  m <- MetricFun(phi = p$phi, theta = p$theta, dtm = d)
  
  m
})

k_m <- do.call(rbind, k_m)

### Vary Lambda ----
lambda_list <- c(50, 100, 250, 500, 1000)
names(lambda_list) <- paste("l", lambda_list, sep="_")

par <- par[[2]] # represents baseline parameters, above

lambda_m <- lapply(lambda_list, function(l){
  d <- SampleDocs(phi = par$phi, theta = par$theta, lambda = l, cpus = 4)
  
  m <- MetricFun(phi = par$phi, theta = par$theta, dtm = d)
  
  m
})

lambda_m <- do.call(rbind, lambda_m)


### Vary V ----
v_list <- c(2000, 5000, 10000, 25000)
names(v_list) <- paste("v", v_list, sep="_")

par <- parallel::mclapply(v_list, 
                          function(v) SimulatePar(D = D, V = v, K = K),
                          mc.cores = parallel::detectCores() - 1)

v_m <- lapply(par, function(p){
  
  d <- SampleDocs(phi = p$phi, theta = p$theta, lambda = lambda)
  
  m <- MetricFun(phi = p$phi, theta = p$theta, dtm = d)
  
  m
})

v_m <- do.call(rbind, v_m)



### Vary D ----
d_list <- c(1000, 2000, 5000, 10000)
names(d_list) <- paste("d", d_list, sep="_")

par <- parallel::mclapply(d_list, 
                          function(d) SimulatePar(D = d, V = V, K = K),
                          mc.cores = parallel::detectCores() - 1)

d_m <- lapply(par, function(p){
  
  d <- SampleDocs(phi = p$phi, theta = p$theta, lambda = lambda)
  
  m <- MetricFun(phi = p$phi, theta = p$theta, dtm = d)
  
  m
})

d_m <- do.call(rbind, d_m)

save(d_m, v_m, k_m, lambda_m, file = "data-derived/r-squared/simulation_metrics.RData")

### NIH DTM ----

nih <- read.csv("data-raw/RePORTER_PRJABS_C_FY2014.csv", stringsAsFactors=FALSE)

s <- sample(seq_len(nrow(nih)), 1000)

nih <- nih$ABSTRACT_TEXT[s]

nih <- stringr::str_conv(nih, "UTF-8")

dtm <- CreateDtm(doc_vec = nih,
                 doc_names = 1:1000,
                 ngram_window = c(1, 1),
                 stopword_vec = stopwords::stopwords("en"),
                 verbose = FALSE)

tf <- TermDocFreq(dtm)

dtm <- dtm[, tf$term[tf$term_freq > 2]]

save(dtm, file = "data-derived/r-squared/nih_sample_dtm.RData")

### TRAIN NIH MDELS ----
k_list <- seq(from = 50, to=500, by=50)

nih_out <- parallel::mclapply(k_list, function(k){
  
  m <- FitLdaModel(dtm = dtm, k = k,
                   iterations = 300,
                   burnin = 250,
                   alpha = 0.1,
                   beta = 0.01, 
                   optimize_alpha = FALSE,
                   calc_likelihood = FALSE,
                   calc_coherence = FALSE,
                   calc_r2 = TRUE)
  
  m$ll <- CalcLikelihood(dtm = dtm, phi = m$phi, theta = m$theta)
  
  data.frame(r2 = m$r2, ll = m$ll, stringsAsFactors = FALSE)
}, mc.cores = parallel::detectCores() - 1)

nih_out <- do.call(rbind, nih_out)

save(nih_out, file = "data-derived/r-squared/nih_out.RData")

### SIM DTM ----
sim_dtm <- SampleDocs(phi = par$d_1000$phi,
                      theta = par$d_1000$theta, 
                      lambda = lambda)

save(sim_dtm, file = "data-derived/r-squared/sim_dtm.RData")

### TRAIN SIM MODELS ----

set.seed(8675309)

k_list <- seq(from = 35, to=70, by=5)

sim_out <- parallel::mclapply(k_list, function(k){
  
  m <- FitLdaModel(dtm = dtm, k = k,
                   iterations = 300,
                   burnin = 250,
                   alpha = 0.1,
                   beta = 0.01, 
                   optimize_alpha = FALSE,
                   calc_likelihood = FALSE,
                   calc_coherence = FALSE,
                   calc_r2 = TRUE)
  
  m$ll <- CalcLikelihood(dtm = dtm, phi = m$phi, theta = m$theta)
  
  data.frame(r2 = m$r2, ll = m$ll, stringsAsFactors = FALSE)
}, mc.cores = parallel::detectCores() - 1)

sim_out <- do.call(rbind, sim_out)

save(sim_out, file = "data-derived/r-squared/sim_out.RData")





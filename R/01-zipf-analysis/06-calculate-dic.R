# This script calculates DIC for models and adds them to the 
# t-statistics tables


### load libraries ---
library(tidyverse)
library(tidylda)
library(Matrix)

### load data ----
dtm_metrics <- read_rds("data-derived/zipf-analysis/dtm-metrics.rds")

model_metrics <- read_rds("data-derived/zipf-analysis/model-metrics.rds")

t_tests <- read_rds("data-derived/zipf-analysis/model-t-tests.rds")

### Augment metrics table with perplexity ----
model_metrics <-
  model_metrics |>
  mutate(
    perplexity = -1 / (Nv * doc_length) * ll
  )



### define a function to calculate DIC ----
calc_dic <- function(m, d, samples = 50){
  
  # sample 50 docs out of posterior
  p <- posterior(m)
  
  theta <- generate(
    x = p,
    matrix = "theta", 
    which = 1:nrow(m$theta),
    times = samples
  ) 
  
  # calc likelihood for `samples` combos
  ll <- by( # for each sample
    theta, 
    INDICES = theta$sample, 
    function(x){
      
      # reshape sample's theta
      t2 <- 
        x |>
        select(-sample) |>
        pivot_wider(
          id_cols = "document", 
          names_from = "topic", 
          values_from = "theta"
        ) |>
        select(-document) |>
        as.matrix()
      
      # get and reshape sample's beta
      b2 <- 
        generate(
          x = p,
          matrix = "beta",
          which = 1:nrow(m$beta),
          times = 1
        ) |> 
        select(-sample) |>
        pivot_wider(
          id_cols = "topic",
          names_from = "token",
          values_from = "beta"
        ) |>
        select(-topic) |>
        as.matrix()
      
      # calculate log likelihood
      ll <- textmineR::CalcLikelihood(
        dtm = d,
        phi = b2,
        theta = t2,
        cpus = 1
      )
      
    }
  ) |>
    unlist() |>
    as.numeric()
  
  # dic is -2 * mean(ll) + 2 * var(ll)
  
  dic <- -2 * mean(ll) + 2 * var(ll)
  
  dic
  
}

### Iterate over models to get each DIC ----
model_names <-
  list.files("data-derived/zipf-analysis/models")

# create batches so that parallel::mclapply doesn't lose results.
num_batches <- parallel::detectCores() - 1

batch_size <- ceiling(length(model_names) / num_batches)

batches <- 
  vector(mode = "list", length = num_batches)

for (j in seq_along(batches)) {
  batches[[j]] <- seq(
    (j - 1) * batch_size + 1,
    min((j - 1) * batch_size + batch_size, length(model_names))
  )
}

batches <- 
  batches |>
  map(function(x) model_names[x])

augmentation <-
  batches |>
  parallel::mclapply(
    function(batch) {
      
      result <-
        batch |>
        map(
          function(filename) {
            # create output table
            out_tab <- 
              str_split(filename, "-")[[1]] 
            
            out_tab <- 
              tibble(
                Nv = as.numeric(out_tab[1]),
                Nd = as.numeric(out_tab[2]),
                Nk = as.numeric(out_tab[3]),
                doc_length = as.numeric(out_tab[4]),
                beta_sum = as.numeric(out_tab[5]),
                alpha_sum = as.numeric(out_tab[6] |> str_replace("\\.rds", "")),
              )
            
            # load model
            models <- read_rds(
              paste0("data-derived/zipf-analysis/models/", filename)
            )
            
            # load dtm
            dat <- read_rds(
              paste0("data-derived/zipf-analysis/simulated-data/", filename)
            )
            
            # calculate relevant stats
            out <-
              models |>
              map(
                function(m){
                  
                  dic <- calc_dic(m = m, d = dat$dtm[, colnames(m$beta)])
                  
                }
              )
            
            # combine with out_tab
            out_tab <- 
              rbind(
                out_tab,
                out_tab,
                out_tab,
                out_tab,
                out_tab
              ) |>
              mutate(
                est_k = 
                  models |> map(function(x) nrow(x$beta)) |> unlist(),
                est_alpha_sum = 
                  models |> map(function(x) nrow(x$beta) * x$alpha) |> unlist() ,
                est_eta_sum = 
                  models |> map(function(x) {
                    if (length(x$eta) > 1) {
                      sum(x$eta)
                    } else {
                      x$eta * ncol(x$beta)
                    }
                  }) |> unlist(),
                est_flat_eta = 
                  models |> map(function(x) {
                    if (length(x$eta) > 1) {
                      FALSE
                    } else {
                      TRUE
                    }
                  }) |> unlist(),
                dic = unlist(out)
              )
            
            # exit
            out_tab
            
          }
        )
      
      # return result
      result |> 
        bind_rows()
      
    }, mc.cores = parallel::detectCores() - 1
  ) |> 
  bind_rows() |>
  mutate( # this ID may not be valid if all cores don't return results
    id = (1:(n() / 5)) |>
      rep(5) |>
      sort()
  )

write_rds(
  augmentation,
  file = "data-derived/zipf-analysis/dic.rds"
)


### Add perplexity to augmentation table ----


augmentation <- 
  augmentation |>
  select(-id) |> # drop potentially invalid ID
  right_join(
    model_metrics |>
      select(
        Nv,
        Nd,
        Nk,
        doc_length,
        beta_sum,
        alpha_sum,
        est_k,
        est_alpha_sum,
        est_eta_sum,
        est_flat_eta,
        perplexity
      )
  ) |> # create a valid ID column
  arrange(
    Nv,
    Nd, 
    Nk,
    doc_length, 
    beta_sum, 
    alpha_sum
  ) |>
  mutate(
    id = (1:(n() / 5)) |>
      rep(5) |>
      sort()
  )


# add indicator if model is correctly specified
augmentation <-
  augmentation |>
  mutate(
    correct_model = 
      Nk == est_k & 
      beta_sum == est_eta_sum & 
      alpha_sum == est_alpha_sum &! 
      est_flat_eta
  )


### Create tables of pairs based on correct/incorrect spec by hyperparam ----

# group by model ID
matched_model_metrics <- 
  augmentation |>
  by(INDICES = augmentation$id, function(x) x)



# If K is too big
nk_big <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_k != x$Nk
      
      if (x$est_k[x$correct_model] < x$est_k[indicator]) {
        out <- tibble(
          id = x$id[1],
          dic_correct = x$dic[x$correct_model],
          dic_est = x$dic[indicator],
          perplexity_correct = x$perplexity[x$correct_model],
          perplexity_est = x$perplexity[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()

# if K is too small
nk_small <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_k != x$Nk
      
      if (x$est_k[x$correct_model] > x$est_k[indicator]) {
        out <- tibble(
          id = x$id[1],
          dic_correct = x$dic[x$correct_model],
          dic_est = x$dic[indicator],
          perplexity_correct = x$perplexity[x$correct_model],
          perplexity_est = x$perplexity[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()


# if eta_sum is too big
eta_sum_big <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_eta_sum != x$beta_sum
      
      if (x$est_eta_sum[x$correct_model] < x$est_eta_sum[indicator]) {
        out <- tibble(
          id = x$id[1],
          dic_correct = x$dic[x$correct_model],
          dic_est = x$dic[indicator],
          perplexity_correct = x$perplexity[x$correct_model],
          perplexity_est = x$perplexity[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()


# if eta_sum is too small
eta_sum_small <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_eta_sum != x$beta_sum
      
      if (x$est_eta_sum[x$correct_model] > x$est_eta_sum[indicator]) {
        out <- tibble(
          id = x$id[1],
          dic_correct = x$dic[x$correct_model],
          dic_est = x$dic[indicator],
          perplexity_correct = x$perplexity[x$correct_model],
          perplexity_est = x$perplexity[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()

# if eta is flat
eta_flat <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_flat_eta
      
      if (TRUE) {
        out <- tibble(
          id = x$id[1],
          dic_correct = x$dic[x$correct_model],
          dic_est = x$dic[indicator],
          perplexity_correct = x$perplexity[x$correct_model],
          perplexity_est = x$perplexity[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()

# if alpha_sum is too big
alpha_sum_big <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_alpha_sum != x$alpha_sum
      
      if (x$est_alpha_sum[x$correct_model] < x$est_alpha_sum[indicator]) {
        out <- tibble(
          id = x$id[1],
          dic_correct = x$dic[x$correct_model],
          dic_est = x$dic[indicator],
          perplexity_correct = x$perplexity[x$correct_model],
          perplexity_est = x$perplexity[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()

# if alpha_sum is too small
alpha_sum_small <-
  matched_model_metrics |>
  map(
    function(x){
      
      # flag for the type of misspecification we want
      indicator <- x$est_alpha_sum != x$alpha_sum
      
      if (x$est_alpha_sum[x$correct_model] > x$est_alpha_sum[indicator]) {
        out <- tibble(
          id = x$id[1],
          dic_correct = x$dic[x$correct_model],
          dic_est = x$dic[indicator],
          perplexity_correct = x$perplexity[x$correct_model],
          perplexity_est = x$perplexity[indicator]
        )
        
        return(out)
        
      } else {
        
        return(NULL)
        
      }
    }
  ) |>
  bind_rows()


### Get t-test on perplexity and DIC for each permutation----

# t-test helper function
t_helper <- function(correct, est) {
  
  my_t <- t.test(est, correct, paired = TRUE)
  
  tibble(
    mean_diff = my_t$estimate,
    ci_95_low = my_t$conf.int[1],
    ci_95_high = my_t$conf.int[2],
    t_stat =  my_t$statistic
  )
  
}


df_helper <- function(df, par) {
  
  out <- list(
    dic = t_helper(df$dic_correct, df$dic_est) |>
      mutate(outcome =  "dic"),
    perplexity = t_helper(df$perplexity_correct, df$perplexity_est) |>
      mutate(outcome =  "perplexity")
  ) |>
    bind_rows() |>
    mutate(
      parameter = par
    )
  
  return(out)
}

t_tests2 <-
  list(
    alpha_sum_big = df_helper(alpha_sum_big, "alpha_sum_big"),
    alpha_sum_small = df_helper(alpha_sum_small, "alpha_sum_small"),
    eta_flat = df_helper(eta_flat, "eta_flat"),
    eta_sum_big = df_helper(eta_sum_big, "eta_sum_big"),
    eta_sum_small = df_helper(eta_sum_small, "eta_sum_small"),
    nk_big = df_helper(nk_big, "nk_big"),
    nk_small = df_helper(nk_small, "nk_small")
  ) |>
  bind_rows()

### write out object(s) for analysis ----

write_rds(
  x = t_tests2,
  file = "data-derived/zipf-analysis/model-t-tests-2.rds"
)

t_tests_combined <- rbind(t_tests, t_tests2)

write_rds(
  x = t_tests_combined,
  file = "data-derived/zipf-analysis/model-t-tests-combined.rds"
)




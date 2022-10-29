# this script loads models and other objects from previous model runs and then 
# produces analysis tables and graphics for the writeup

library(tidyverse)
library(tidytext)
library(tidylda)

### Define a very important variable ----
### load relevant data ----
pars <- read_rds("data-derived/tlda-sims/pop-pars.rds")

hdist <- read_rds("data-derived/tlda-sims/hdist.rds")


### Get time series of hdist differences ----

hdiff <- lapply(
  hdist,
  function(a_values) { # loop over different values of a
    parallel::mclapply( # loop over the different data generating distributions
      a_values, 
      function(par_values) {
        
        out <- vector(mode = "list", length = length(par_values))
        
        out[[1]] <- par_values[[1]] / par_values[[1]]
        
        for (j in 2:length(out)) { # loop over each model by period (skip the first)
          out[[j]] <- (par_values[[j]] - par_values[[j - 1]]) / par_values[[1]]
        }
        
        out
      }, mc.cores = parallel::detectCores() - 1
    )
  }
)

write_rds(
  hdiff,
  file = "data-derived/tlda-sims/hdiff.rds"
)

### Reformat hdiff for analysis ----

# make into a tidy tibble
hdiff_clean <- 
  hdiff %>% map(
    function(a_values) { # for each value of a
      a_values %>%
        parallel::mclapply(
          function(par_values) { # for each data generating parameter 
            
            out <- vector(mode = "list", length = length(par_values))
            
            for (j in seq_along(out)) { # for each time period
              out[[j]] <- par_values[[j]] %>%
                as_tibble() %>%
                mutate(actual_topic = rownames(par_values[[j]])) %>%
                pivot_longer(
                  cols = matches("^[0-9]+$"),
                  names_to = "fitted_topic",
                  values_to = "hdiff"
                ) %>%
                mutate(period = j)
            }
            
            out <- out %>% bind_rows()
            
            out
          }, mc.cores = parallel::detectCores() - 1
        )
    }
  )


# add generating parameter values
hdiff_clean <- 
  hdiff_clean %>%
  map(
    function(a_values) {
      for (j in seq_along(a_values)) {
        a_values[[j]] <- 
          a_values[[j]] %>%
          mutate(
            model_num = pars[[j]]$model_num, # corresponds to combo of {doc_length, beta_sum, alpha_sum, alpha_shape}
            doc_length = pars[[j]]$doc_length,
            beta_sum = pars[[j]]$par$beta_sum,
            alpha_sum = pars[[j]]$par$alpha_sum,
            alpha_shape = pars[[j]]$par$alpha_shape
          )
      }
      
      a_values
    }
  )

# collapse into a single tibble
hdiff_clean <- 
  hdiff_clean %>%
  map(
    function(a_values) {
      a_values %>%
        map(
          function(par_values){
            par_values %>% bind_rows()
          }
        ) %>%
        bind_rows()
    }
  ) 

for (j in seq_along(hdiff_clean)) {
  
  hdiff_clean[[j]] <-
    hdiff_clean[[j]] %>%
    mutate(
      a = str_sub(names(hdist)[j], start = 3, end = -1) %>% as.numeric()
    )
  
}

hdiff_clean <- 
  hdiff_clean %>%
  bind_rows()

# save hdiff_clean
write_rds(
  hdiff_clean,
  file = "data-derived/tlda-sims/hdiff_clean.rds"
)

### check for model convergence ----
# Do a test in two stages:
# 1. Do a t-test on the last 20 periods of each series of models for each pair of
#    actual topic and fitted_topic. The topic estimates converged if we reject 
#    the null that mean hdiff is 0 over that period at the 0.05 significance level.
# 2. Do a binomal test on the whole model. Under the null, we'd expect to reject
#    5% of topics by random chance. Model has converged if we reject the null of
#    a 5% "success rate" of rejecting topics. (Corresponds to 1.5 topics.) Use
#    0.05 significance level here too.

# t-test on all topics
get_t_test_p_val <- function(x, mu, alternative = "two.sided") {
  result <- t.test(x, mu = mu, alternative = alternative)
  
  result$p.value
}

t_table <-
  hdiff_clean %>%
  group_by(
    a,
    model_num, 
    doc_length,
    beta_sum,
    alpha_sum,
    alpha_shape,
    fitted_topic,
    actual_topic
  ) %>%
  filter(period > 80) %>%
  summarize(
    t_test_p = get_t_test_p_val(hdiff, mu = 0)
  )

# test of proportions for each model
get_bin_test_p_val <- function(success, fail, p) {
  
  result <- prop.test(
    x = success,
    n = success + fail,
    p = p,
    alternative = "greater"
  )
  
  result$p.value
}

bin_table <-
  t_table %>%
  ungroup() %>%
  mutate(
    topic_converge = t_test_p <= 0.05
  ) %>%
  group_by(
    a,
    model_num, 
    doc_length,
    beta_sum,
    alpha_sum,
    alpha_shape,
  ) %>%
  summarize(
    num_converge = sum(topic_converge),
    num_not_converge = sum(1 - topic_converge)
  ) %>%
  group_by(
    a,
    model_num, 
    doc_length,
    beta_sum,
    alpha_sum,
    alpha_shape,
  ) %>%
  summarize(
    bin_test_p = get_bin_test_p_val(
      success = num_converge,
      fail = num_not_converge,
      p = 0.05
    )
  )

bin_regression <- 
  glm(
    bin_test_p <= 0.05 ~ a + doc_length + beta_sum + alpha_sum, 
    data = bin_table,
    family = binomial("logit")
  )

# save objects for convergence analysis
save(
  t_table,
  bin_table,
  bin_regression,
  file = "data-derived/tlda-sims/convergence-analysis.RData"
)

# clean up to preserve memory space
rm(
  dtms,
  hdiff,
  hdiff_clean
)

gc()


### reformat hdist for analysis ----
hdist_clean <- 
  hdist %>% map(
    function(a_values) { # for each value of a
      a_values %>%
        parallel::mclapply(
          function(par_values) { # for each data generating parameter 
            
            out <- vector(mode = "list", length = length(par_values))
            
            for (j in seq_along(out)) { # for each time period
              out[[j]] <- par_values[[j]] %>%
                as_tibble() %>%
                mutate(actual_topic = rownames(par_values[[j]])) %>%
                pivot_longer(
                  cols = matches("^[0-9]+$"),
                  names_to = "fitted_topic",
                  values_to = "hdist"
                ) %>%
                mutate(period = j)
            }
            
            out <- out %>% bind_rows()
            
            out
          }, mc.cores = parallel::detectCores() - 1
        )
    }
  )

# add generating parameter values
hdist_clean <- 
  hdist_clean %>%
  map(
    function(a_values) {
      for (j in seq_along(a_values)) {
        a_values[[j]] <- 
          a_values[[j]] %>%
          mutate(
            model_num = pars[[j]]$model_num, # corresponds to combo of {doc_length, beta_sum, alpha_sum, alpha_shape}
            doc_length = pars[[j]]$doc_length,
            beta_sum = pars[[j]]$par$beta_sum,
            alpha_sum = pars[[j]]$par$alpha_sum,
            alpha_shape = pars[[j]]$par$alpha_shape
          )
      }
      
      a_values
    }
  )

# collapse into a single tibble
hdist_clean <- 
  hdist_clean %>%
  map(
    function(a_values) {
      a_values %>%
        map(
          function(par_values){
            par_values %>% bind_rows()
          }
        ) %>%
        bind_rows()
    }
  ) 

for (j in seq_along(hdist_clean)) {
  
  hdist_clean[[j]] <-
    hdist_clean[[j]] %>%
    mutate(
      a = str_sub(names(hdist_clean)[j], start = 3, end = -1) %>% as.numeric()
    )
  
}

hdist_clean <- 
  hdist_clean %>%
  bind_rows()

# save hdist_clean
save(
  hdist_clean,
  file = "data-derived/tlda-sims/hdist_clean.RData"
)

### get topic matches ----

# filter by combinations where the model converged
hdist_clean_converged <- 
  bin_table %>%
  ungroup() %>%
  filter(bin_test_p <= 0.05) %>%
  select(
    a,
    model_num
  ) %>%
  inner_join(
    hdist_clean,
    by = c("a" = "a", "model_num" = "model_num")
  )



# find mateches by averaging over last 20 periods for converged topics
match_by_avg_periods <- 
  hdist_clean_converged %>%
  filter(period > 80) %>%
  group_by(
    fitted_topic,
    a,
    model_num,
    doc_length,
    beta_sum,
    alpha_sum,
    alpha_shape
  ) %>%
  summarize(
    best_match = actual_topic[which.min(mean(hdist))],
    hdist = min(mean(hdist))
  )

save(
  match_by_avg_periods,
  file = "data-derived/tlda-sims/topic_matches.RData"
)

### See if average hdiff for matched topics is negative ----
hdiff_clean <- read_rds("data-derived/tlda-sims/hdiff_clean.rds")


# As above, do a test in two stages:
# 1. Do a t-test on the last 98 periods of each series of models for each pair of
#    actual topic and fitted_topic. The topic estimates converged towards the
#    target topic if we reject the null that the mean hdiff over that period is
#    greater than zero at the 0.05 significance level.
# 2. Do a binomal test on the whole model. Under the null, we'd expect to reject
#    5% of topics by random chance. Model has converged towards the target if we 
#    reject the null of a 5% "success rate" of rejecting topics. Use the 0.05
#    significance level here too.

# filter out the first two periods and non-matches
# then do a t-test for each topic match's derivative
matched_hdiff_t_test <-
  hdiff_clean %>%
  filter(period > 2) %>%
  inner_join(
    match_by_avg_periods,
    by = c(
      "a" = "a",
      "model_num" = "model_num",
      "doc_length" = "doc_length",
      "beta_sum" = "beta_sum",
      "alpha_sum" = "alpha_sum",
      "alpha_shape" = "alpha_shape",
      "fitted_topic" = "fitted_topic",
      "actual_topic" = "best_match"
    )
  ) %>%
  group_by(
    a,
    model_num,
    actual_topic,
    fitted_topic,
    doc_length,
    beta_sum,
    alpha_sum,
    alpha_shape
  ) %>%
  summarise(
    mean_hdiff = mean(hdiff),
    t_test_p = get_t_test_p_val(x = hdiff, mu = 0, alternative = "greater")
  )

# get the binomial test
get_bin_exact_p_val <- function(success, fail, p) { # need exact b/c smaller sample sizes
  
  result <- binom.test(
    x = success,
    n = success + fail,
    p = p,
    alternative = "greater"
  )
  
  result$p.value
}

matched_bin_table <-
  matched_hdiff_t_test %>%
  ungroup() %>%
  mutate(
    topic_converge = t_test_p <= 0.05
  ) %>%
  group_by(
    a,
    model_num, 
    doc_length,
    beta_sum,
    alpha_sum,
    alpha_shape,
  ) %>%
  summarize(
    num_converge = sum(topic_converge),
    num_not_converge = sum(1 - topic_converge)
  ) %>%
  group_by(
    a,
    model_num, 
    doc_length,
    beta_sum,
    alpha_sum,
    alpha_shape,
  ) %>%
  summarize(
    bin_test_p = get_bin_exact_p_val(
      success = num_converge,
      fail = num_not_converge,
      p = 0.05
    )
  )


save(
  matched_hdiff_t_test,
  matched_bin_table,
  file = "data-derived/tlda-sims/convergence-direction-analysis.RData"
)

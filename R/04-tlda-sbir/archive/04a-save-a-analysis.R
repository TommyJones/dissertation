# running 04-sbir-a-analysis.R made me run out of memory
# take temporary saved results and do something with them


# load libraries
library(tidyverse)
library(textmineR)
library(tidylda)
library(tidytext)


a_range <- seq(0.4, 1.6, by = 0.2)

year_range <- 1984:2004

# load from file to save time for now
load("data-derived/tlda-sbir/sbir_models-by-award-year.RData")

m_1983 <- sbir_models[[1]]

rm(sbir_models)

gc()


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

model_metrics <- 
  parallel::mclapply(
    a_range,
    function(a) {
      
      metric_list <- vector(mode = "list", length = length(year_range))
      
      names(metric_list) <- year_range
      
      model_path <- paste0("data-derived/tlda-sbir/models-a-", a, "-", year_range[1], ".rds")
      
      m <- read_rds(model_path)
      
      metric_list[[1]] <- tibble(
        a = a,
        year = year_range[1],
        likelihood = mean(m$log_likelihood$log_likelihood[51:200]),
        r2 = m$r2,
        mean_coh = mean(m$summary$coherence, na.rm = TRUE),
        var_coh = var(m$summary$coherence, na.rm = TRUE),
        mean_hellinger = mean(get_hellinger_delta(m_1983, m))
      )
      
      for(j in 2:length(metric_list)) {
        
        model_path <- paste0("data-derived/tlda-sbir/models-a-", a, "-", year_range[j], ".rds") 
        
        m2 <- read_rds(model_path)
        
        metric_list[[j]] <- 
          tibble(
            a = a,
            year = year_range[j],
            likelihood = mean(m2$log_likelihood$log_likelihood[51:200]),
            r2 = m2$r2,
            mean_coh = mean(m2$summary$coherence, na.rm = TRUE),
            var_coh = var(m2$summary$coherence, na.rm = TRUE),
            mean_hellinger = mean(get_hellinger_delta(m, m2))
          )
        
        m <- m2
      }
      
      metric_list <-
        metric_list |> bind_rows()
      
      metric_list

    }, mc.cores = parallel::detectCores() - 1
  ) |>
  bind_rows()

save(
  model_metrics,
  file = "data-derived/tlda-sbir/sbir-a-metrics.RData"
)


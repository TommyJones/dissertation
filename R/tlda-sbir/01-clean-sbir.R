# script cleans raw SBIR data for later use

# load libraries
library(tidyverse)
library(textmineR)
library(tidylda)
library(lubridate)

# load raw data
sbir <- read_csv("https://data.www.sbir.gov/awarddatapublic/award_data.csv")

colnames(sbir) <- 
  colnames(sbir) %>%
  tolower() %>%
  str_replace_all(" +", "_")

sbir <- 
  sbir %>%
  mutate(
    sbir_id = 1:nrow(sbir),
    proposal_award_date = as_date(proposal_award_date, format = "%m/%d/%Y"),
    contract_end_date = as_date(contract_end_date, format = "%m/%d/%Y")
  ) 


# pull out text columns
sbir_text <- sbir %>%
  select(
    sbir_id,
    award_title,
    abstract
  ) %>%
  mutate(
    award_title = str_conv(award_title, "UTF-8"),
    abstract = str_conv(abstract, "UTF-8")
  )

save(
  sbir,
  sbir_text,
  file = "data-derived/tlda-sbir/sbir-clean.RData"
)
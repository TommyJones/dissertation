# This script makes some sample data for developing a sampler in Rust
# use the NIH sample data from textmineR and then downsample number of 
# documents and tokens to have a small data set that can be printed on
# screen if needed for inspection

library(tidyverse)
library(Matrix)

d <- textmineR::nih_sample_dtm

dim(d)

d2 <- d[1:10, ]

d2 <- d2[, colSums(d2) > 3]

dim(d2)

triplet <- as(d2, "dgTMatrix") 
  
# note i and j in triplet are 0 indexed (as in Rust)
triplet <- tibble(i = triplet@i, j = triplet@j, v = triplet@x)

vocab <- colnames(d2)

doc_names <- rownames(d2)

# get the raw documents as well. Maybe play with hugging face tokenizers later
documents <- textmineR::nih_sample[1:10, c("APPLICATION_ID", "ABSTRACT_TEXT")]

# write data out
if (! dir.exists("data_derived/rust-sample-data"))
  dir.create("data_derived/rust-sample-data")

if (! dir.exists("data_derived/rust-sample-data/documents"))
  dir.create("data_derived/rust-sample-data/documents")

write_csv(
  x = triplet,
  file = "data_derived/rust-sample-data/triplet-dtm.csv",
  na = "",
  append = FALSE,
  col_names = FALSE,
  quote_escape = "none"
  )

write_csv(
  x = tibble(x= vocab),
  file = "data_derived/rust-sample-data/vocab.txt",
  na = "",
  append = FALSE,
  col_names = FALSE,
  quote_escape = "none"
)

write_csv(
  x = tibble(x = doc_names),
  file = "data_derived/rust-sample-data/doc-names.txt",
  na = "",
  append = FALSE,
  col_names = FALSE,
  quote_escape = "none"
)

for (j in 1:nrow(documents))
  write_file(
    x = documents[j, "ABSTRACT_TEXT"],
    file = paste0("data_derived/rust-sample-data/documents/", 
                  documents[j, "APPLICATION_ID"],
                  ".txt")
  )



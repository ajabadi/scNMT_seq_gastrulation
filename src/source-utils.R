library(here)
sapply( list.files(here('src/utils'), pattern = "*.R", full.names = TRUE), source)

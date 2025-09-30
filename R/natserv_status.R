library(natserv)
library(tidyverse)
library(dplyr)
library(DBI)

# Search for species by common name or scientific name
res <- ns_search_spp(text = "Abies balsamea", page = 0, per_page = 5)
best_result <- res$results[1,]

nations <- best_result$nations[[1]]
ca <- nations %>% 
  filter(nationCode == 'CA')
subnations <- ca[, "subnations"][[1]]
qc <- subnations %>% 
  filter(subnationCode == 'QC')

exotic = qc[,'exotic']
native = qc[,'native']

View(qc)
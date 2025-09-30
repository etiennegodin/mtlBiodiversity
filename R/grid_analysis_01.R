library(tidyverse)
library(dplyr)
library(DBI)
library(duckdb)
library(arrow)
library(sf) ## for geom



con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/db/mtlbio.duckdb", read_only = TRUE)
dbListTables(con)
df <- dbReadTable(con, "grid_sjoin")   # read into R

df <- df %>%
  filter(!is.na(grid_id))

#Remove ids from user sampling bias filter
loaded_ids <- readRDS("R/ids_samplingBias.rds")
loaded_ids <- as.numeric(loaded_ids)

filtered_df <- df[!df$gbifID %in% loaded_ids, ]
print(nrow(df) - nrow(filtered_df))
print('Removed')





#Shannon index helper function
shannon <- function(x) {
  freq <- table(x)          #count of each species
  p <- freq / sum(freq)     # proportions
  -sum(p*log(p))            # H
}

simpson <- function(x) {
  freq <- table(x)
  p <- freq / sum(freq)
  1 - sum(p^2)                    # Simpson diversity
}

df_eco <- filtered_df %>% 
  group_by(grid_id) %>% 
  summarise(
    observations = n(),
    species_richness = n_distinct(species),
    shannon_index = shannon(species),
    simpson_index = simpson(species),
    genus_richness = n_distinct(genus),
    family_richness = n_distinct(family),
    order_richness = n_distinct(order),
    class_richness = n_distinct(class),
    phylum_richness = n_distinct(phylum),
    n_users = n_distinct(recordedBy),
    H_user = shannon(recordedBy),
    D_user = simpson(recordedBy)
    
    
  )

view(df_eco)

write_parquet(df_eco, "data/processed/grid.parquet")



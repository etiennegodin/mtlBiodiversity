library(tidyverse)
library(dplyr)
library(DBI)
library(duckdb)
library(arrow)
library(sf) ## for geom

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

filter_df <- function(df){
  
  df <- df %>%
    filter(!is.na(grid_id))
  
  #Remove ids from user sampling bias filter
  loaded_ids <- readRDS("R/ids_samplingBias.rds")
  loaded_ids <- as.numeric(loaded_ids)
  
  filtered_df <- df[!df$gbifID %in% loaded_ids, ]
  print(nrow(df) - nrow(filtered_df))
  print('Removed')
  return(filtered_df)
  
}

analysis <- function(df, group_col){
  
  df_filtered <- filter_df(df)

  df_analysis <- df_filtered %>% 
    group_by({{group_col}}) %>% 
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
      D_user = simpson(recordedBy),
      rpo = species_richness/observations,
      user_diveristy = species_richness/n_users
  
    )
  print('Analysis done')
  return(df_analysis)
}

con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/db/mtlbio.duckdb", read_only = TRUE)
dbListTables(con)

df_grid <- dbReadTable(con, "grid_sjoin")   # read into R
df_quartiers <- dbReadTable(con, "quartiers_sjoin")   # read into R
df_parks <- dbReadTable(con, "parks_sjoin")   # read into R

grid_analysis_df <- analysis(df_grid, grid_id)
quartiers_analysis_df <- analysis(df_quartiers, qrt_id)
parks_analysis_df <- analysis(df_parks, park_id)

write_parquet(grid_analysis_df, "data/processed/grid.parquet")
write_parquet(quartiers_analysis_df, "data/processed/quartiers.parquet")
write_parquet(parks_analysis_df, "data/processed/parks.parquet")



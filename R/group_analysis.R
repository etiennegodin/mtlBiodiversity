if (requireNamespace("renv", quietly = TRUE)) {
  renv::activate()
}

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

filter_df <- function(df, group_col){

  # Remove rows where not sjoined in grouping col 
  df <- df %>%
    filter(!is.na(group_col))
  
  #Drop cols where rows are na on more than 60%
  threshold <- 0.6
  na_fracs <- colMeans(is.na(df))
  df <- df[,na_fracs <= threshold]

  #Remove ids from user sampling bias filter
  loaded_ids <- readRDS("R/ids_samplingBias.rds")
  loaded_ids <- as.numeric(loaded_ids)
  
  filtered_df <- df[!df$gbifID %in% loaded_ids, ]
  print(nrow(df) - nrow(filtered_df))
  print('Removed')
  return(filtered_df)
  
}

analysis <- function(df, group_col){
  
  df_filtered <- filter_df(df, group_col)

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
#snakemake objects 
db_file <- snakemake@params[['db_name']]
table <- snakemake@params[['table']]
group_col <- snakemake@params[['group_table']]

con <- dbConnect(duckdb::duckdb(), dbdir = db_file, read_only = TRUE)
df <- dbReadTable(con, table)   # read into R
grid_analysis_df <- analysis(table, group_col)
write_parquet(grid_analysis_df, "snakemake@output[['out]])




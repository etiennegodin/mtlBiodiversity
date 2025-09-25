library(tidyverse)
library(DBI)
library(duckdb)
library(arrow)
library(sf) ## for geom

con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/db/mtlbio.duckdb", read_only = TRUE)
dbListTables(con)
df <- dbReadTable(con, "quartiers_sjoin")   # read into R

df <- df %>%
  filter(!is.na(qrt_id))
view(df)


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

df_eco <- df %>% 
  group_by(qrt_id) %>% 
  summarise(
    observations = n(),
    species_richness = n_distinct(species),
    shannon_index = shannon(species),
    simpson_index = simpson(species),
    genus_richness = n_distinct(genus),
    family_richness = n_distinct(family),
    order_richness = n_distinct(order),
    class_richness = n_distinct(class),
    phylum_richness = n_distinct(phylum)
  )

view(df_eco)

write_parquet(df_eco, "data/processed/quartiers.parquet")



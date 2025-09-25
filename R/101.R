library(DBI)
library(duckdb)
library(tidyverse)

con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/db/mtlbio.duckdb", read_only = TRUE)
dbListTables(con)
df <- dbReadTable(con, "quartiers_sjoin")   # read into R

view(df)
glimpse((df))

?unique

df %>% 
  group_by(qrt_id) %>% 
  summarise(
    species_richness = n_distinct(species),
    genus_richness = n_distinct(genus),
    family_richness = n_distinct(family),
    order_richness = n_distinct(order),
    class_richness = n_distinct(class)
    
    
  )


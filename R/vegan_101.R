library(vegan)
library(dplyr)
library(tidyr)
library(duckdb)

filter_df <- function(df, group_id){
  
  df <- df %>%
    filter(!is.na(.data[[group_id]]))
  
  #Remove ids from user sampling bias filter
  loaded_ids <- readRDS("R/ids_samplingBias.rds")
  loaded_ids <- as.numeric(loaded_ids)
  
  filtered_df <- df[!df$gbifID %in% loaded_ids, ]
  print(nrow(df) - nrow(filtered_df))
  print('Removed')
  return(filtered_df)
  
}

# Data prep 
con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/db/mtlbio.duckdb", read_only = TRUE)
dbListTables(con)

table <- "quartiers_sjoin"
group_id <- 'qrt_id'

df <- dbReadTable(con, table)   # read into R

df_filtered <- filter_df(df, group_id)

df_counts <- df_filtered %>% 
  add_count(.data[[group_id]], name = 'n_obs') %>% 
  filter(n_obs >= 100 )

comm <- df_counts %>% 
  count(.data[[group_id]], species) %>% 
  pivot_wider(names_from = species, values_from = n, values_fill = 0)


comm <- as.data.frame(comm)
rownames(comm) <- comm[[group_id]]
comm <- comm[, -1]  # drop park_value column

parks_obsrvations <- rowSums(comm)
parks_obsrvations


## Alpha diversity ##

# Shannon index
diversity(comm, index = "shannon")

# Simpson index
diversity(comm, index = "simpson")

# Richness (number of species)
specnumber(comm)

## Sampling effort correction ##
raremax <- min(rowSums(comm))  # smallest sample size

# Plot rarefaction curves up to the smallest sample size
rarecurve(comm, step = round(raremax/20), sample = raremax,
          xlab = "Number of individuals", ylab = "Expected species",
          label = TRUE)

# Compute rarefied richness at that common sample size (numeric)
rarefied_richness <- rarefy(comm, sample = raremax)
park_rich_rare <- data.frame(site = rownames(comm),
           sample_size = rowSums(comm),
           rarefied_richness = rarefied_richness)

## BETA Diversity ## 

# Bray-Curtis dissimilarity
dist_comm <- vegdist(rarefied_richness, method = "bray")
#View(as.matrix(dist_comm))

# Hierarchical clustering
clust <- hclust(dist_comm)
plot(clust)


# Bray-Curtis NMDS
set.seed(123)  # NMDS is random
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100)
# Basic plot
plot(nmds)

# Add site labels
ordiplot(nmds, type = "n")
orditorp(nmds, display = "sites", col = "blue", cex = 1.2)
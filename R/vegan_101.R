library(vegan)
library(dplyr)
library(tidyr)
library(duckdb)
library(arrow)
filter_df <- function(df, group_id){
  
  df <- df %>%
    filter(!is.na(.data[[group_id]]))
  
  #Remove ids from user sampling bias filter
  loaded_ids <- readRDS("R/ids_samplingBias.rds")
  loaded_ids <- as.numeric(loaded_ids)
  
  filtered_df <- df[!df$gbifID %in% loaded_ids, ]
  
  filtered_df <- filtered_df %>%
    filter(phylum == 'Tracheophyta')
  
  print(nrow(df) - nrow(filtered_df))
  print('Removed')
  return(filtered_df)
  
}

# Data prep 
con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/db/mtlbio.duckdb", read_only = TRUE)
dbListTables(con)

table <- "parks_sjoin"
group_id <- 'park_id'
cluster_count = 3
df <- dbReadTable(con, table)   # read into R

?dbDisconnect(con, shutdown = TRUE)
df_filtered <- filter_df(df, group_id)



df_counts <- df_filtered %>% 
  add_count(.data[[group_id]], name = 'n_obs') %>% 
  filter(n_obs >= 10 ) %>% 
  count(.data[[group_id]], species)

group_ids <- unique(df_counts[[group_id]])
df_out <- data.frame(park_id = group_ids )

comm <- df_counts %>% 
  pivot_wider(names_from = species, values_from = n, values_fill = 0)

comm <- as.data.frame(comm)
rownames(comm) <- comm[[group_id]]
comm <- comm[, -1]  # drop park_value column


## Alpha diversity ##

# Shannon index
df_out$shannon_index <- unlist(diversity(comm, index = "shannon"))

# Simpson index
df_out$simpson_index <- unlist(diversity(comm, index = "simpson"))

# Richness (number of species)
df_out$species_richness <- unlist(specnumber(comm))

n_obs = rowSums(comm)
df_out$n_obs <- unlist(n_obs)

## Sampling effort correction ##
raremax <- min(rowSums(comm))  # smallest sample size

# Plot rarefaction curves up to the smallest sample size
rarecurve(comm, step = round(raremax/5), sample = raremax,
          xlab = "Number of individuals", ylab = "Expected species",
          label = TRUE)

# Compute rarefied richness at that common sample size (numeric)
rarefied_richness <- rarefy(comm, sample = raremax)
df_out$rarefied_richness <- unlist(rarefied_richness)


## BETA Diversity ## 

# Bray-Curtis dissimilarity
dist_comm <- vegdist(rarefied_richness, method = "bray")
#View(as.matrix(dist_comm))

# Hierarchical clustering
clust <- hclust(dist_comm, method = 'average')
# Cut tree into k clusters (e.g., 3)
clusters_labels <- cutree(clust, k = cluster_count)
# Turn into dataframe

df_out$clusters <- unlist(clusters_labels)

k_values <- 2:6
avg_sil <- numeric(length(k_values))

library(cluster)  # for silhouette

for (i in seq_along(k_values)) {
  k <- k_values[i]
  clusters <- cutree(clust, k = k)
  sil <- silhouette(clusters, dist(dist_comm))
  avg_sil[i] <- mean(sil[, 3])
}

# Plot average silhouette vs k
plot(k_values, avg_sil, type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Average silhouette width",
     main = "Silhouette-based k selection")



library(ggplot2)
library(reshape2)

# Convert distance matrix to dataframe
bc_mat <- as.matrix(dist_comm)
bc_df <- melt(bc_mat, varnames = c("site1","site2"), value.name = "bray")

plot(clust, labels = rownames(comm), main = "Cluster Dendrogram (Bray-Curtis)")
rect.hclust(clust, k = cluster_count, border = "red")  # highlight clusters

# Bray-Curtis NMDS
# NMDS with Bray–Curtis
set.seed(123)
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100)

# Add NMDS coordinates to dataframe
nmds_points <- as.data.frame(scores(nmds, display = "sites"))
nmds_points$park_id <- group_ids

df_out <- df_out %>% 
  left_join(nmds_points, by = group_id)

# Plot NMDS with cluster colors
ggplot(df_out, aes(x = NMDS1, y = NMDS2, color = factor(clusters))) +
  geom_point(size=4) +
  theme_minimal() +
  labs(color="Cluster", title="NMDS of Parks (Bray-Curtis)") +
  stat_ellipse(level=0.95)


adonis_res <- adonis2(dist_comm ~ clusters_labels + n_obs, data = df_out)



bd <- betadisper(dist_comm, df_out$cluster)
anova(bd)  # traditional ANOVA on distances to centroids
permutest(bd)  # permutation test (more robust)

# Boxplot of distances to centroid
boxplot(bd, main="Multivariate Dispersion (beta diversity)", ylab="Distance to centroid")

# Or plot in NMDS space
plot(bd)


write_parquet(df_out, "data/processed/vegan.parquet")

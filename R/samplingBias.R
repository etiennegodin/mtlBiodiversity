library(dbscan)
library(geosphere)  # for distance matrix
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(duckdb)
library(arrow)
library(sf)
library(DBI)



con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/db/mtlbio.duckdb", read_only = TRUE)
dbListTables(con)
df <- dbReadTable(con, "gbif_raw")   # read into R

#View(df)

df_top <- df %>% 
  group_by(identifiedBy) %>% 
  summarise(obs_count = n()) %>% 
  arrange(desc(obs_count))
View(df_top)

users <- filter(df_top, obs_count > 1000 )
users <- users[!is.na(users$identifiedBy),]
users <- users$identifiedBy
#user <- filter(users, )
View(users)


cluster_stats <- function(df, coords_utm)
{
  # ---------------------------------------------------------
  # 3. Summarise cluster stats
  # ---------------------------------------------------------
  df$x <- coords_utm[,1]
  df$y <- coords_utm[,2]
  
  cluster_stats <- df %>%
    group_by(cluster) %>%
    summarise(
      n_obs = n(),
      n_species = n_distinct(species),
      n_users = n_distinct(identifiedBy),
      # spread = mean distance to centroid
      spread_mean_m = {
        centroid <- colMeans(cbind(x, y))
        mean(sqrt((x - centroid[1])^2 + (y - centroid[2])^2))
      },
      spread_max_m = {
        centroid <- colMeans(cbind(x, y))
        max(sqrt((x - centroid[1])^2 + (y - centroid[2])^2))
      }
    )
  
  #print(cluster_stats)
  flags <- cluster_stats %>%
    mutate(
      flag_home = (n_obs > 100 & n_users == 1 & spread_mean_m < 30 ),
      flag_chaining = (n_obs > 100 & n_users == 1 & spread_mean_m < 30 & spread_max_m > 200 )
    )
  
  selected_clusters <- c()
  for (i in 1:nrow(flags))
  {
    if (flags$flag_home[i] == "TRUE")
    {
      selected_clusters <- c(selected_clusters, flags$cluster[i] )
    }
  }
  
  return(selected_clusters)

}

db_scan_func <- function(df)
{
  
  coords <- df[, c('decimalLongitude', 'decimalLatitude')]
  df_sf <- st_as_sf(df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  df_utm <- st_transform(df_sf, 32618)  # pick your UTM zone
  coords_utm <- st_coordinates(df_utm)
  
  db <- dbscan(coords_utm, eps = 50, minPts = 2)
  #db_h <- hdbscan(coords_utm, minPts = 5)
  df$cluster <- db$cluster

  clusters <- cluster_stats(df, coords_utm)
  #print(clusters)

  if (length(clusters) > 0 )
  {
    print(paste(length(clusters),'Clusters flagged'))
    for (i in 1:nrow(df))
    {
      if (df$cluster[i] %in% clusters){
        ids_to_remove <<- c(ids_to_remove, df$gbifID[i])
        #print(ids_to_remove)
      
      }
    }
  }
  #print(length(ids_to_remove))
  
}


# Initialize empty vector to store ids
ids_to_remove <- c()

for (u in users)
{
  print(u)
  df_temp <- filter(df, identifiedBy == u)

  ids <-db_scan_func(df_temp)
}
print(length(ids_to_remove))

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




for (u in users)
{
  print(paste("Hello", u))
  next
  df_temp <- filter(df, identifiedBy == u)
  #?distm
  coords <- df_temp[, c('decimalLongitude', 'decimalLatitude')]
  df_sf <- st_as_sf(df_temp, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  df_utm <- st_transform(df_sf, 32618)  # pick your UTM zone
  coords_utm <- st_coordinates(df_utm)

  #db <- dbscan(coords_utm, eps = 50, minPts = 5)
  db_h <- hdbscan(coords_utm, minPts = 5)
  df_temp$cluster <- db_h$cluster
}



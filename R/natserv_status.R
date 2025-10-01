library(natserv)
library(tidyverse)
library(dplyr)
library(DBI)
library(duckdb)


get_status <- function(species_name){
  #print(species_name)
  exotic <- NA
  native <- NA
  b_result <- NULL
  # Search for species by common name or scientific name
  res <- ns_search_spp(text = species_name, page = 0, per_page = 5)
  #Check if results are provided
  if (length(res) > 1) {
    results <- res$results
  } else {
    stop("No results found")
  }
  
  if (length(results) == 0) stop("Results is empty")

  #Check if results 
  #Iterate over each species to find the one matching with provided name 
  for (r in 1:nrow(results)) {
    if (results$scientificName[r] == species_name) {
      b_result <- results %>% filter(scientificName == species_name)
      break
    }
  }
    
  #If no matching name, return NA as unsure 
  if (is.null(b_result) || nrow(b_result) == 0) stop("No matching species")
  
  # If right species name, extract exotic/native values
  nations <- b_result$nations[[1]]
  if (nrow(nations) == 0) stop("No national data")
    
  ca <- nations %>% filter(nationCode == 'CA')
  if (nrow(ca) == 0) stop("No data for Canada")

  
  subnations <- ca[, "subnations"][[1]]
  qc <- subnations %>% filter(subnationCode == 'QC')
  
  if (nrow(qc) != 0) {
    exotic <- qc$exotic
    native <- qc$native
  } else {
    exotic <- ca$exotic
    native <- ca$native
  }

  # Assign status
  if (exotic != native) {
    if (native == TRUE) return(0)
    if (exotic == TRUE) return(1)
  } else {
    return(2)
  }
} 

safe_apply <- function(species_list) {
  pbapply::pblapply(seq_along(species_list), function(i) {
    s <- species_list[[i]]
    tryCatch({
      list(
        index = i,
        species = s,
        status = get_status(s),
        error = NA_character_
      )
    }, error = function(e) {
      list(
        index = i,
        species = s,
        status = NA,
        error = e$message
      )
    })
  }) %>%
    lapply(as.data.frame) %>%
    do.call(rbind, .)
}

con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/db/mtlbio.duckdb", read_only = TRUE)
dbListTables(con)
df <- dbReadTable(con, "grid_sjoin")   # read into R

df <- df %>%
  filter(!is.na(grid_id))

df <- df %>%
  filter(!is.na(species))

species_list <- unique(df$species)
#View(species_list)
#species_list <- species_list[1:300]
#species_list <-c("Polygonia satyrus", "Caligo telamonius","Myscelus assaricus")

df_status <- safe_apply(species_list)
View(df_status)

write_parquet(df_status, "data/processed/speciesQcStatus.parquet")
write.csv(df_status, file = "data/processed/speciesQcStatus.csv", row.names = FALSE)


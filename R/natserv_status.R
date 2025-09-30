library(natserv)
library(tidyverse)
library(dplyr)
library(DBI)
library(duckdb)


get_status <- function(species_name){
  exotic <- NA
  native <- NA
  b_result <- NA
  res <- NA
  # Search for species by common name or scientific name
  res <- ns_search_spp(text = species_name, page = 0, per_page = 5)
  if (length(res)> 1)
  {
    results <- res$results
    #View(results)
  }
  for (r in 1:nrow(results))
  {
    if (results$scientificName[r] != species_name){
      next
    }
    else if (results$scientificName[r] == species_name){
      b_result <- results %>% 
        filter(scientificName == species_name)
      #View(b_result)
    }
  if (!length(b_result)> 1){
    return(NA)
  }
  else
  {
    nations <- b_result$nations[[1]]
    if (nrow(nations) == 0){
      return(NA)
    }
    
    ca <- nations %>% 
      filter(nationCode == 'CA')
    print(ca)
    
    #View(ca)
    subnations <- ca[, "subnations"][[1]]
    #View(subnations)
    qc <- subnations %>% 
      filter(subnationCode == 'QC')
    if(!nrow(qc) == 0){
      exotic = qc[,'exotic']
      native = qc[,'native']
    }
    else{
      exotic = ca[,'exotic']
      native = ca[,'native']
    }

    if (is.na(exotic) & is.na(native)){
      return(NA)
    }
    else{
      if (exotic != native)
      {
        if (native == TRUE)
        {
          return(0)
        }
        else if ( exotic == TRUE)
        {
          return(1)
        }
      }
      else
      {
        return(2)
      }
    }

    
  }

  }
  
}

con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/db/mtlbio.duckdb", read_only = TRUE)
dbListTables(con)
df <- dbReadTable(con, "grid_sjoin")   # read into R

df <- df %>%
  filter(!is.na(grid_id))

species_list <- unique(df$species)
View(species_list)

#species_list <- species_list[1:20]

species_list <-c("Rumex sanguineus")


statuses <- c()
for (s in species_list){
  print(s)
  status = get_status(s)
  statuses <<- c(statuses, status)
  
}
View(statuses)

df_status <- data.frame(
  species = species_list,
  status = statuses
)
View(df_status)
write.csv(df_status, file = "data/speciesQcStatus.csv", row.names = FALSE)


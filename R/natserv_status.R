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
  #Check if results are provided
  if (length(res)> 1)
  {
    #Extract first list as df results 
    results <- res$results
  }
  #Check if results list is empty
  if (length(results) == 0 ){
    stop("Results is empty")   # will trigger an error
  }

  #Check if results 
  #Iterate over each species to find the one matching with provided name 
  for (r in 1:nrow(results))
  {
    if (results$scientificName[r] != species_name){
      next
      #If not species_name, skip to next iteration 
    }
    else if (results$scientificName[r] == species_name){
      #If item in list matches, store as best results
      b_result <- results %>% 
        filter(scientificName == species_name)
      #View(b_result)
    }
    
  #If no matching name, return NA as unsure 
  if (!length(b_result)> 1){
    stop('No species data match to provided species')   # will trigger an error
  }
  # If right species name, extract exotic/native values
  else
  {
    nations <- b_result$nations[[1]]
    if (nrow(nations) == 0){
      stop('No national data')   # will trigger an error
    }
    
    ca <- nations %>% 
      filter(nationCode == 'CA')
    if (nrow(ca) == 0){
      stop('No data for canada')   # will trigger an error
    }
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
      stop('No exotic or native data')   # will trigger an error
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

df <- df %>%
  filter(!is.na(species))

species_list <- unique(df$species)
View(species_list)

#species_list <- species_list[1:2]

#species_list <-c("Tiarella cordifolia")


statuses <- c()
for (s in species_list){
  print(s)
  
  status <- tryCatch({
    get_status(s)
  }, error = function(e) {
    message("Error for element ", s, ": ", e$message)
    return(NA)  # or NA, or just skip
  })

  statuses <<- c(statuses, status)
  
}
View(statuses)
View(results)

df_status <- data.frame(
  species = species_list,
  status = statuses
)
View(df_status)

write_parquet(df_status, "data/speciesQcStatus.parquet")
write.csv(df_status, file = "data/speciesQcStatus.csv", row.names = FALSE)


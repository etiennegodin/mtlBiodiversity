#etl

# Folder structure
- Gather data and dump in raw/data
	- data/raw/parks 
	- data/raw/gbif
	- data/raw/city
- Make data/raw read only by python 

### Only data/raw has variable file names - otherwise set programmatically in interim 


# Park


| parkd_id | park_name | park_type1 | park_type1 |
| -------- | --------- | ---------- | ---------- |
|          |           |            |            |

##  uniformise_park_data 
1. Iterate over `*.shp` in data/raw/parks
2. Check if file interim is done, --force 
3. For each, check if data fields exists 
	1. If not prompt user with `park_x_.shp` columns 
	2. Loop over each table field to select from `park_x.shp` columns 
	3. Save as a config 
		1. `park_x_.shp` OBJECTID -> parkd_id 


## clean_park_data 

1. 
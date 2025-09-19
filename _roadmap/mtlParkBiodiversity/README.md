
# 🍄🐦🌳 Montreal Biodiversity Dashboard🌳🐦🍄

Explore the richness of Montreal’s urban nature through an interactive dashboard that visualizes biodiversity across the city’s parks and neighborhoods. This project demonstrates geospatial data analysis, ecological metrics computation, and interactive visualization in Python.

## Project Highlights

- **Data Integration:** Combined species occurrence data from the Global Biodiversity Information Facility (GBIF) with Montreal's open data.
    
- **Geospatial Processing:** Cleaned and prepared shapefiles from Montreal's datasets. Then performed spatial joins to assign species observations to Create grid over the city's extent t.
	
- **Biodiversity Metrics:** Calculated key metrics such as species richness, abundance, and diversity indices for parks and neighborhoods.
	
- **Time series:** Calculated key metrics such as species richness, abundance, and diversity indices for parks and neighborhoods.
	
- **Interactive Visualization:** Built a Python dashboard with maps, charts, and filters to explore biodiversity patterns across Montreal.

## Impact & Insights

- Identify parks and neighborhoods with the highest species diversity.
	
- Enable the general public to discover Montreal's hidden diversity and richness
	
- Explore patterns of urban biodiversity across Montreal.
	
- Support research, city planning, and conservation efforts with accessible, location-based biodiversity data.
	
- Demonstrate the integration of open biodiversity data with geospatial analysis for actionable insights.

## How It Works

1. Choose your explore mode either 
	1. Your personal address
	2. Choose your neighborhood from Montreal's list
	3. Choose a park to explore 
2. Run mbio.exe with command  `full` or either more granularly:
	1. `geo` - run park and city shapefile preprocess
	2. `gbif` - run spatial join 
	3. `metric` - process metrics 
	4. `app` - launch dashboard on local server 
    


### Sources:
GBIF.org (05 September 2025) GBIF Occurrence Download  https://doi.org/10.15468/dl.7jaynu
https://donnees.montreal.ca/en/dataset/quartiers
https://donnees.montreal.ca/en/dataset/grands-parcs-parcs-d-arrondissements-et-espaces-publics
CREATE OR REPLACE TABLE parks AS

SELECT park_name,
COUNT(DISTINCT gbifID) AS observation, -- dsda-- 
COUNT(DISTINCT kingdom) AS kingdom_richness,
COUNT(DISTINCT family) AS family_richness,
COUNT(DISTINCT genus) AS genus_richness,
COUNT(DISTINCT species) AS species_richness,


ANY_VALUE(park_geom) AS park_geom,
    

FROM data d

WHERE park_name IS NOT NULL
GROUP BY park_name
ORDER BY species_richness DESC

;

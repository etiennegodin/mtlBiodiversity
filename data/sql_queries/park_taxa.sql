CREATE OR REPLACE TABLE park_taxa AS



SELECT park_name,
kingdom, class, genus, species,

COUNT (species) AS species_count,
ANY_VALUE(park_id) AS park_geom,
ANY_VALUE(park_geom) AS park_geom,


FROM data
GROUP BY park_name, kingdom, class, genus, species,
ORDER BY park_name, kingdom, class, genus, species
;

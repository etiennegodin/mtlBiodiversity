CREATE OR REPLACE TABLE park_taxa AS



SELECT park_name,kingdom, class, genus, species,
COUNT (species) AS species_count,
FROM data
GROUP BY park_name, kingdom, class, genus, species,
ORDER BY park_name, kingdom, class, genus, species
;

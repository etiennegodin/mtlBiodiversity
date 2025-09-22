CREATE OR REPLACE TABLE {{new_table_name}} AS

SELECT {{group_by}},
COUNT(DISTINCT gbifID) AS observation, -- dsda-- 
COUNT(DISTINCT kingdom) AS kingdom_richness,
COUNT(DISTINCT family) AS family_richness,
COUNT(DISTINCT genus) AS genus_richness,
COUNT(DISTINCT species) AS species_richness,

ANY_VALUE({{group_by_id}}),
ANY_VALUE(geom),
    

FROM {{table_name}}

WHERE {{{group_by}}} IS NOT NULL
GROUP BY {{group_by}}
ORDER BY species_richness DESC

;

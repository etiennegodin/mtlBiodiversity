
CREATE OR REPLACE VIEW park_shannon_index AS
WITH species_count AS
(SELECT park_name, species, 
COUNT(*) AS species_count,
FROM data
WHERE park_name IS NOT NULL
GROUP BY park_name, species,
ORDER BY park_name
),

park_totals AS(

SELECT park_name,
SUM(species_count) AS total_count
FROM species_count
GROUP BY park_name,
),

proportions AS(

SELECT s.park_name,
        s.species,
        s.species_count,
        t.total_count,
        s.species_count * 1.0 / ( t.total_count *1.0 )AS p
        
FROM species_count s
JOIN park_totals t
ON s.park_name = t.park_name
)


SELECT park_name,

    -SUM(p*LN(p)) AS shannon_index
FROM proportions
GROUP BY park_name;
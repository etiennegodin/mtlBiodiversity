CREATE OR REPLACE TABLE taxa_groups AS

WITH kingdom_richness AS(
    SELECT park_name, kingdom,
    COUNT(*) AS kingdom_count
    FROM data
    GROUP BY park_name, kingdom,

),

class_richness AS(
    SELECT park_name, class,
    COUNT(*) AS class_count
    FROM data
    GROUP BY park_name, class,
),

taxa_groups AS(
        SELECT k.park_name,
        k.kingdom,
        k.kingdom_count,
        c.class_count,
        c.class

    FROM kingdom_richness k
    JOIN class_richness c
    ON k.park_name = c.park_name

)

SELECT *
FROM taxa_groups
WHERE park_name IS NOT NULL

;

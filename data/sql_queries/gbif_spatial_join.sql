                CREATE OR REPLACE TABLE gbif_with_parks AS
                SELECT g.gbifID,
                g.occurrenceID,
                g.kingdom,
                g.phylum,
                g.class,
                g.order,
                g.family,
                g.genus,
                g.species,
                g.taxonRank,
                g.scientificName,
                g.eventDate,
                g.day,
                g.month,
                g.year,
                g.taxonKey,
                g.identifiedBy,
                g.basisOfRecord,
                g.license,
                g.recordedBy,
                g.geom,

                {park_fields_sql}

                p.geom AS park_geom,
                FROM gbif g
                LEFT JOIN parks p
                    ON g.maxx >= ST_XMIN(p.geom)
                    AND g.minx <= ST_XMAX(p.geom)
                    AND g.maxy >= ST_YMIN(p.geom)
                    AND g.miny <= ST_YMAX(p.geom)
                    AND ST_Within(g.geom, p.geom) -- Spatial join predicate
                ;
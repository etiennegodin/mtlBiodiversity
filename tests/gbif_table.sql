CREATE OR REPLACE TABLE observations AS
                SELECT f.gbifID,
                f.occurrenceID,
                f.kingdom,
                f.phylum,
                f.class,
                f.order,
                f.family,
                f.genus,
                f.species,
                f.taxonRank,
                f.scientificName,
                f.eventDate,
                f.day,
                f.month,
                f.year,
                f.taxonKey,
                f.basisOfRecord,
                f.license,
                f.recordedBy,
                
                ST_Point(decimalLongitude, decimalLatitude) AS geom,
                FROM 'data\interim\gbif\_test_gbif_data.parquet' AS f
                LIMIT 10 
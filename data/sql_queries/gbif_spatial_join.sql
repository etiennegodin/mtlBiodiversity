CREATE OR REPLACE TABLE {{output_table_name}}  AS 
                SELECT 
                
                {{left_fields}} 
                {{right_fields}} 
                r.geom AS right_geom,

                FROM {{left_table_name}} l 
                LEFT JOIN {{right_table_name}} r 
                    ON l.maxx >= r.minx
                    AND l.minx <= r.maxx
                    AND l.maxy >= r.miny
                    AND l.miny <= r.maxy
                    AND ST_Intersects(l.geom, r.geom) -- Spatial join predicate
                ;

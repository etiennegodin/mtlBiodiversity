CREATE OR REPLACE TABLE {{output_table_name}}  AS 
                SELECT 
                
                {{left_fields}} 
                {{right_fields}} 
                r.geom AS right_geom,

                FROM {{left_table_name}} l 
                LEFT JOIN {{right_table_name}} r 
                    ON l.maxx >= ST_XMIN(r.geom)
                    AND l.minx <= ST_XMAX(r.geom)
                    AND l.maxy >= ST_YMIN(r.geom)
                    AND l.miny <= ST_YMAX(r.geom)
                    AND ST_Within(l.geom, r.geom) -- Spatial join predicate
                ;

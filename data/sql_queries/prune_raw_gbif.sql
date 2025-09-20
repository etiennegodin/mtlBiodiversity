ALTER TABLE {{output_table_name}} DROP COLUMN geom;
ALTER TABLE {{output_table_name}} DROP COLUMN minx;
ALTER TABLE {{output_table_name}} DROP COLUMN miny;
ALTER TABLE {{output_table_name}} DROP COLUMN maxx;
ALTER TABLE {{output_table_name}} DROP COLUMN maxy;


ALTER TABLE {{output_table_name}} RENAME COLUMN right_geom TO geom;
ALTER TABLE {{output_table_name}} RENAME COLUMN minx_1 TO minx;
ALTER TABLE {{output_table_name}} RENAME COLUMN miny_1 TO miny;
ALTER TABLE {{output_table_name}} RENAME COLUMN maxx_1 TO maxx;
ALTER TABLE {{output_table_name}} RENAME COLUMN maxy_1 TO maxy;



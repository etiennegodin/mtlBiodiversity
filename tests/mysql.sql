ALTER TABLE _test_occurences_grid DROP COLUMN geom;
ALTER TABLE _test_occurences_grid DROP COLUMN minx;
ALTER TABLE _test_occurences_grid DROP COLUMN miny;
ALTER TABLE _test_occurences_grid DROP COLUMN maxx;
ALTER TABLE _test_occurences_grid DROP COLUMN maxy;


ALTER TABLE _test_occurences_grid RENAME COLUMN right_geom TO geom;
ALTER TABLE _test_occurences_grid RENAME COLUMN minx_1 TO minx;
ALTER TABLE _test_occurences_grid RENAME COLUMN miny_1 TO miny;
ALTER TABLE _test_occurences_grid RENAME COLUMN maxx_1 TO maxx;
ALTER TABLE _test_occurences_grid RENAME COLUMN maxy_1 TO maxy;


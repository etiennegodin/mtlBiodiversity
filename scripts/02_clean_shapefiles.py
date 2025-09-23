from mtlBio.dataprep.clean_shp import main
print('#'*25, '02_clean_shapefiles', '#'*25,'\n' )

main(input_file = snakemake.input[0], output_file = snakemake.output[0], crs = snakemake.params.crs)

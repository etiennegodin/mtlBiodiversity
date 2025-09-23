
from mtlBio.dataprep.gbif import convert_gbif_csv

convert_gbif_csv(snakemake.input[0], snakemake.output[0], limit= snakemake.params.limit)
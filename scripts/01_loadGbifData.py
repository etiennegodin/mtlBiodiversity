
from mtlBio.dataprep.gbif import import_gbif_csv

import_gbif_csv(snakemake.input[0], snakemake.output[0], limit= snakemake.params.limit)
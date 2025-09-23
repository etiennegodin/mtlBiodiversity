
from mtlBio.dataprep.gbif import import_gbif_csv

print('#'*25, '01_loadGbifData', '#'*25, '\n' )

import_gbif_csv(snakemake.input[0], snakemake.output[0], limit= snakemake.params.limit)
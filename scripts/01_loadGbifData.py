
from mtlBio.dataprep.importGbif import main

print('#'*25, '01_loadGbifData', '#'*25, '\n' )

main(snakemake.input[0], snakemake.output[0], limit= snakemake.params.limit)
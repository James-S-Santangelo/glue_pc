rsync -vuar -P ../../../results/ \
    --exclude '*.zarr/*'  \
    santang3@nia-datamover1.scinet.utoronto.ca:/scratch/n/nessrobe/santang3/toronto_gwsd/results

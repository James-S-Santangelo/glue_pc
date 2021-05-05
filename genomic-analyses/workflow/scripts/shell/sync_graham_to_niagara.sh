rsync -vuar -P \
    --exclude 'qc' \
    --exclude 'figures' \
    --exclude 'merged' \
    --exclude 'paired' \
    --exclude 'unpaired' \
    --exclude 'final/s_*' \
    ../../../results \
    santang3@nia-datamover1.scinet.utoronto.ca:/scratch/n/nessrobe/santang3/glue/glue-paper1/genomic-analyses

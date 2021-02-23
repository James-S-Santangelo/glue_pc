# ANGSD
rsync -vuar -P \
    --prune-empty-dirs \
    --include='*/' \
    --include='*.sfs' \
    --exclude='*' \
    santang3@graham.computecanada.ca:/scratch/santang3/glue-low1/snakemake/results/single_sample_sfs/ \
    ../../data/single_sample_sfs


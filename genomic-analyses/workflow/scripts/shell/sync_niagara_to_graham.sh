rsync -vuar -P \
    --dry-run \
    santang3@nia-datamover1.scinet.utoronto.ca:/scratch/n/nessrobe/santang3/glue/glue-paper1/genomic-analyses/workflow/logs \
    ../../logs/

rsync -vuar -P \
    --dry-run \
    santang3@nia-datamover1.scinet.utoronto.ca:/scratch/n/nessrobe/santang3/glue/glue-paper1/genomic-analyses/workflow/slurm_logs \
    ../../slurm_logs/

rsync -vuar -P \
    --dry-run \
    santang3@nia-datamover1.scinet.utoronto.ca:/scratch/n/nessrobe/santang3/glue/glue-paper1/genomic-analyses/results/angsd \
    ../../../results/angsd/

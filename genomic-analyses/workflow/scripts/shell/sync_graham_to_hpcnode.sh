rsync -vuar -P \
    santang3@gra-dtn1.computecanada.ca:/home/santang3/scratch/glue/glue-paper1/genomic-analyses/workflow/logs/ \
    ../../logs

rsync -vuar -P \
    santang3@gra-dtn1.computecanada.ca:/home/santang3/scratch/glue/glue-paper1/genomic-analyses/workflow/slurm_logs/ \
    ../../../slurm_logs

rsync -vuar -P \
    santang3@gra-dtn1.computecanada.ca:/home/santang3/scratch/glue/glue-paper1/genomic-analyses/results/ \
    ../../../results

This directory contains text files used in the pipeline:

- [chromosome_file.txt](./chromosome_file.txt): Simple text file with names of 16 chromosomes in white clover reference genome.
- [low1_libraryConcentrations](./low1_libraryConcentrations.csv): File used to get sample IDs to parse raw read files. __note__: Can probably delete this and call the file from the [sequencing-prep](../../sequencing-prep/) directly
- [toronto_samples_for_glue](./toronto_samples_for_glue.txt): Plants from Toronto to include in pipeline and downsample to 1X. __Note__: Should probably incorporate raw reads for these samples into the main pipeline instead of calling the BAMs for these files. Can then delete this and incorporate into the main sample sheet. 
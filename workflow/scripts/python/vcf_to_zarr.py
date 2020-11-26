# Script to convert SNP VCF to Zarr database

import sys

import numcodecs
import allel

def vcf_to_zarr(vcf_in, outpath):
    """Convert on-disk VCF to on-disk Zarr database using
    scikit-allele and zarr modules

    Zarr database written to same directory as input VCF
    
    Args:
        vcf_in (str): Path to input VCF on disk
        outpath (str): Path to derectory to output Zarr db
        chrom (str): Chromosome for which Zarr database should be created

    Returns:
        None
    """
    # Rename 'numalt' field. Required by Zarr to distinguish `NUMALT` from `numalt`
    # `numalt` is automatically computed by scikit-allel
    rename_dict = {'variants/numalt':'variants/numalt_sci'}

    # Use vcf_to_zarr function from scikit-allel to create zarr database
    # Currently optimized for biallelic SNP VCF but easy to extend functionality
    allel.vcf_to_zarr(
            input=vcf_in,
            output=outpath,
            overwrite=True,
            rename_fields=rename_dict,
            fields='*',
            alt_number=1,
            compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False)
            )

with open(snakemake.log[0], 'w') as f:
    sys.stderr = sys.stdout = f
    vcf_in = snakemake.input[0]
    outpath = snakemake.output[0]
    vcf_to_zarr(vcf_in, outpath)

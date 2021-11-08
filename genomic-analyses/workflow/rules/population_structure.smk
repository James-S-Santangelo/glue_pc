# Population structure analyses. 

################
#### GLOBAL ####
################

rule pcangsd_allSamples:
    """
    Perform PCA using genome-wide 4fold dengenerate sites using all samples from all cities.
    """
    input:
        rules.concat_angsd_gl.output
    output:
        '{0}/pcangsd/allSamples/allSamples_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
    log: 'logs/pcangsd_allSamples/allSamples_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 10
    params:
        out = '{0}/pcangsd/allSamples/allSamples_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 12000,
        time = '03:00:00'
    shell:
        """
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            -iter 5000 \
            > {log}
        """

###################
#### ADMIXTURE ####
###################

rule pcangsd_byCity:
    """
    Perform PCA by city using genome-wide 4fold dengenerate sites and estimate admixture proportions
    """
    input:
        rules.angsd_gl_byCity.output.gls
    output:
        cov = '{0}/pcangsd/by_city/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
        pi = '{0}/pcangsd/by_city/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.pi.npy'.format(POP_STRUC_DIR)
    log: 'logs/pcangsd_byCity/{city}_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 6
    params:
        out = '{0}/pcangsd/by_city/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '02:00:00'
    shell:
        """
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            -admix \
            -admix_seed 42 \
            -iter 5000 \
            > {log}
        """
                
rule pop_structure_done:
    """
    Generate empty flag file signaling successful completion of PCAngsd
    """
    input:
        expand(rules.pcangsd_allSamples.output, site = '4fold', maf = ['0.05']),
        expand(rules.pcangsd_byCity.output, site = '4fold', maf = '0.05', city = CITIES)
    output:
        '{0}/population_structure.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """

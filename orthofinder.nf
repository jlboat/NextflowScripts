#!/usr/bin/env nextflow

params.directory = "$baseDir/data/references"

process orthoFinder {
    publishDir "results"

    clusterOptions = { "-l select=1:ncpus=12:mem=5gb,walltime=48:00:00" }
    module 'singularity'

    script:
    """
    singularity run -B /scratch2 ~/singularity_containers/orthofinder.simg -f ${params.directory} -t 12 -a 1 \
-S diamond -M msa -A mafft -T fasttree -o ${baseDir}/results/
    """
}

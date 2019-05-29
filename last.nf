#!/usr/bin/env nextflow

/*
 * Defines pipeline parameters in order to specify the 
 * read pairs by using the command line options
 */
params.chromosomes = "$baseDir/data/references/Sbicolor.Chr*.softmasked.fa"
params.genome_m = "$baseDir/data/references/Zmays_493_APGv4.softmasked.fa"
   
genome_file_m = file(params.genome_m)
sobic_chromosomes = Channel
                    .fromPath(params.chromosomes)
                    .map { file -> tuple(file.baseName, file) }

index_chromosomes = Channel
                    .fromPath(params.chromosomes)
                    .map { file -> tuple(file.baseName, file) }

/*
 * Index chromosomes
 */
process indexChr {
    publishDir "$baseDir/data/references"

    input:
        set dataset_id, file(chromosome) from index_chromosomes

    output:
        file "${dataset_id}.???" into indices

    clusterOptions = { "-l select=1:ncpus=2:mem=12gb,walltime=05:00:00" }
    module 'singularity'

    script:
    """
    singularity run -B /scratch2,/zfs ~/singularity_containers/last.simg lastdb -c -P 2 ${dataset_id} ${chromosome}
    """
}

/*
 * Last alignment
 */
process alignment {
    publishDir "results"

    input:
        set dataset_id, file(chromosome) from sobic_chromosomes
        file "${dataset_id}.???" from indices
        file genome_file_m
 
    output:
      file "Zmays_onto_Sbicolor.Chr*.softmasked.maf"

    clusterOptions = { "-l select=1:ncpus=4:mem=12gb,walltime=05:00:00" } 
    module 'singularity'
 
    script:
    """
    singularity run -B /scratch2,/zfs ~/singularity_containers/last.simg lastal -P4 $baseDir/data/references/${dataset_id} ${genome_file_m} > Zmays_onto_${dataset_id}.maf
    """
} 


#!/usr/bin/env nextflow

/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/fastq/*fastq.gz"
   
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .map { file -> tuple(file.simpleName, file) }
    .into { reads; fastqc_reads }

/*
 * QC raw reads
 */
process preFastqc {
    publishDir "$baseDir/results/fastqc/pre/"

    input:
        set val(dataset_id), file(r1) from fastqc_reads

    output:
        file "*_fastqc.{zip,html}" into raw_fastqc_results

    module 'singularity'

    script:
    """
singularity run -B /zfs,/scratch2 ~/singularity_containers/fastqc.simg -q ${r1}
    """
}

/*
 * Clean the reads
 */
process cleanReads {
    publishDir "$baseDir/data/fastq/trimmed/"

    input:
      set val(dataset_id), file(forward) from reads
 
    output:
      set dataset_id, file("${dataset_id}.TRIMMED.fq.gz") into trimmed_reads
 
    module 'singularity'
    clusterOptions = '-l select=1:ncpus=4:mem=5gb,walltime=08:00:00'
 
    script:
    """
singularity run -B /zfs,/scratch2 ~/singularity_containers/fastp.simg --length_required=80 --thread=4 --in1=${forward} --out1=${dataset_id}.TRIMMED.fq.gz
    """
} 

/*
 * QC trimmed reads
 */
process postFastqc {
    publishDir "results/fastqc/post/"

    input:
        set val(dataset_id), file(r1) from trimmed_reads

    output:
        file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    module 'singularity'

    script:
    """
singularity run -B /zfs,/scratch2 ~/singularity_containers/fastqc.simg -q ${r1} 
    """
}




#!/usr/bin/env nextflow

/*
 * Parameters for execution
 */
params.genome = "$baseDir/Zmays_493_APGv4.fa"
params.gtf = "$baseDir/Zmays_493_RefGen_V4.gene_exons.gff3"
   
/*
 * The reference genome file
 */
genome_gtf = file(params.gtf)  
genome_fna = file(params.genome)


/*
 * Index Genome using STAR
 */
process indexReference {
    publishDir "$baseDir/data/references/"

    input:
        file genome_gtf
        file genome_fna

    clusterOptions = {"-l select=1:ncpus=6:mem=90gb,walltime=25:00:00"}
    module 'singularity'

    script:
    """
singularity run -B /scratch2 ~/singularity_containers/STAR.simg STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${baseDir} \
--genomeFastaFiles ${genome_fna} \
--sjdbGTFfile ${genome_gtf} \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 100
    """

}


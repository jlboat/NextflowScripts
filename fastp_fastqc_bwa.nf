#!/usr/bin/env nextflow

/*
 * Defines pipeline parameters in order to specify the 
 * read pairs by using the command line options
 */
params.reads = "$baseDir/data/*{1,2}.fastq.gz"
params.genome = "$baseDir/data/hg38.fa"
   
/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { fastqc_read_pairs }

genome_file = file(params.genome)

/*
 * Step 1. Clean the reads
 */
process cleanReads {
    publishDir "results"

    input:
      set dataset_id, file(forward), file(reverse) from read_pairs
 
    output:
      set dataset_id, file("${dataset_id}1.TRIMMED.PAIRED.fq.gz"), file("${dataset_id}2.TRIMMED.PAIRED.fq.gz") into fastqc_paired_fastq
      set dataset_id, file("${dataset_id}1.TRIMMED.UNPAIRED.fq.gz"), file("${dataset_id}2.TRIMMED.UNPAIRED.fq.gz") into unpaired_fastq
      set dataset_id, file("${dataset_id}1.TRIMMED.PAIRED.fq.gz"), file("${dataset_id}2.TRIMMED.PAIRED.fq.gz") into bwa_paired_fastq

    clusterOptions = { "-l select=1:ncpus=1:mem=3gb,walltime=04:00:00" } 
    module 'singularity'
 
    script:
    """
    singularity run -B /scratch2 ~/singularity_containers/fastp.simg \
-i $forward -I $reverse \
-o ${dataset_id}1.TRIMMED.PAIRED.fq.gz -O ${dataset_id}2.TRIMMED.PAIRED.fq.gz \
--unpaired1 ${dataset_id}1.TRIMMED.UNPAIRED.fq.gz --unpaired2 ${dataset_id}2.TRIMMED.UNPAIRED.fq.gz \
--thread ${task.cpus} --length_required 40 --detect_adapter_for_pe \
-j ${dataset_id}.json -h ${dataset_id}.html 
    """
} 

/*
 * Step 2. QC trimmed reads
 */
process postFastqc {
    publishDir "results"

    input:
        set dataset_id, file(r1), file(r2) from fastqc_paired_fastq

    output:
        file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    clusterOptions = { "-l select=1:ncpus=1:mem=2gb,walltime=04:00:00" }
    module 'singularity'

    script:
    """
    singularity run -B /scratch2 ~/singularity_containers/fastqc.simg -q $r1 $r2
    """
}

/*
 * Step 3. BWA MEM alignment
 */
process bwa {
    publishDir "results"

    input:
        set dataset_id, file(r1), file(r2) from bwa_paired_fastq
        file genome from genome_file

    output:
        file "*.bam" into bam_files

    clusterOptions = { "-l select=1:ncpus=2:mem=6gb,walltime=04:00:00" }
    module 'singularity'
    module 'anaconda'

    script:
    """
# load samtools from aligners env
source activate aligners
singularity run -B /scratch2 ~/singularity_containers/bwa.simg mem -t ${task.cpus} ${params.genome} ${r1} ${r2} | samtools view -b - > ${dataset_id}.bam
    """

}


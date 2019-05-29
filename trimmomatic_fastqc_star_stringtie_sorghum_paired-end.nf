#!/usr/bin/env nextflow

/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/reads/*fq.gz"
params.genome = "$baseDir/data/references/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
params.gtf = "$baseDir/data/references/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gtf"
   
/*
 * The reference genome file
 */
genome_gtf = file(params.gtf)  
genome_fna = file(params.genome)

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { fastqc_read_pairs }

/*
 * Step 1. Clean the reads
 */
process cleanReads {
    publishDir "$baseDir/data/reads/"

    input:
      set dataset_id, file(forward), file(reverse) from read_pairs
 
    output:
      set dataset_id, file("${dataset_id}_R1.TRIMMED.PAIRED.fq.gz"), file("${dataset_id}_R2.TRIMMED.PAIRED.fq.gz") into fastqc_paired_fastq
      set dataset_id, file("${dataset_id}_R1.TRIMMED.UNPAIRED.fq.gz"), file("${dataset_id}_R2.TRIMMED.UNPAIRED.fq.gz") into unpaired_fastq
 
    clusterOptions = { "--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=3gb --cpus-per-task=1" }
    module 'trimmomatic'
 
    script:
    """
    trimmomatic \
    PE \
    -trimlog /dev/null \
    -threads ${task.cpus} \
    $forward $reverse \
    ${dataset_id}_R1.TRIMMED.PAIRED.fq.gz ${dataset_id}_R1.TRIMMED.UNPAIRED.fq.gz \
    ${dataset_id}_R2.TRIMMED.PAIRED.fq.gz ${dataset_id}_R2.TRIMMED.UNPAIRED.fq.gz \
    ILLUMINACLIP:/apps/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10:2:true \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MAXINFO:40:0.3 \
    MINLEN:40
    """
} 

/*
 * QC trimmed reads
 */
process postFastqc {
    publishDir "results"

    input:
    set dataset_id, file(r1), file(r2) from fastqc_paired_fastq

    output:
    file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    clusterOptions = { "--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=2gb --cpus-per-task=1" }
    module 'fastqc'

    script:
    """
    fastqc -q $r1 $r2
    """
}


/*
 * Index Genome using STAR
 */
process indexReference {
    publishDir "$baseDir/data/references/"

    input:
        file genome_gtf
        file genome_fna

    clusterOptions = {--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=2gb --cpus-per-task=2}
    module 'star'

    script:
    """
    STAR --runThreadN 2
         --runMode genomeGenerate
         --genomeDir ${baseDir}/data/references/
         --genomeFastaFiles ${genome_fna}
         --sjdbGTFfile ${genome_gtf}
         --sjdbOverhang 100
    """

}

/*
 * Align reads using STAR
 */
process alignReads {
    publishDir "$baseDir/data/alignments/"

    input:
        set dataset_id, file(r1), file(r2) from fastqc_paired_fastq

    output:
        set dataset_id, file(bam_file) from bam_files

    clusterOptions = {--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=2gb --cpus-per-task=2}
    module 'star'

    script:
    """
    STAR --runThreadN 2
         --genomeDir ${baseDir}/data/references/
         --genomeFastaFiles ${genome_fna}
         --readFilesIn $r1 $r2
    """

}

/*
 * Assembles the transcript using "stringtie"
 * and publish the transcript output files into the `results` folder
 */
process makeTranscript {
     
    input:
      set dataset_id, file(bam_file) from bam_files
    
    output:
      set '*' into transcripts
  
    clusterOptions = { "--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=3gb --cpus-per-task=1 --no-requeue" }
    module 'stringtie'

    """
    stringtie -e -B -p ${task.cpus} -j 2 -o ${dataset_id}_transcripts.gtf -G ${genome_gtf} ${bam_file}

#    mkdir ${baseDir}/results/${dataset_id}
#    mv ${baseDir}/results/*.ctab ${baseDir}/results/${dataset_id}/
#    mv ${baseDir}/results/*.gtf ${baseDir}/results/${dataset_id}/
    """
}


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
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .map { file -> tuple(file.simpleName, file) }
    .into { reads; fastqc_reads }

/*
 * QC raw reads
 */
process preFastqc {
    publishDir "results/fastqc/pre/"

    input:
    set val(dataset_id), file(r1) from fastqc_reads

    output:
    file "*_fastqc.{zip,html}" into raw_fastqc_results

    clusterOptions = { "--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=2gb --cpus-per-task=1" }
    module 'fastqc'

    script:
    """
    fastqc -q ${r1}
    """
}

/*
 * Clean the reads
 */
process cleanReads {
    publishDir "$baseDir/results/reads/"

    input:
      set val(dataset_id), file(forward) from reads
 
    output:
      set dataset_id, file("${dataset_id}.TRIMMED.fq.gz") into trimmed_reads
      file("${dataset_id}.TRIMMED.fq.gz") into trimmed_reads_2
 
    clusterOptions = { "--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=3gb --cpus-per-task=1" }
    module 'trimmomatic'
 
    script:
    """
    trimmomatic \
SE \
-trimlog /dev/null \
-threads ${task.cpus} \
${forward} ${dataset_id}.TRIMMED.fq.gz \
ILLUMINACLIP:/apps/trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:30:10:2:true \
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
    publishDir "results/fastqc/post/"

    input:
    set val(dataset_id), file(r1) from trimmed_reads

    output:
    file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    clusterOptions = { "--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=2gb --cpus-per-task=1" }
    module 'fastqc'

    script:
    """
    fastqc -q ${r1} 
    """
}


/*
 * Index Genome using STAR
 */
process indexReference {
    publishDir "$baseDir/data/references/"

    input:
        file genome_gtf from genome_gtf
        file genome_fna from genome_fna

    output:
        file '*' into indexed

    clusterOptions = { "--account=icbrbi --qos=icbrbi --time=5:00:00 --mem-per-cpu=5gb --cpus-per-task=4" }
    module 'star'

    script:
    """
    STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir ${baseDir}/data/references/ \
--genomeFastaFiles ${genome_fna} \
--sjdbGTFfile ${genome_gtf} \
--sjdbOverhang 100
    """

}

/*
 * Align reads using STAR
 */
process alignReads {
    publishDir "$baseDir/results/alignments/"

    input:
        file '*' from indexed
        file r1 from trimmed_reads_2

    output:
        file "*.bam" into bam_files

    clusterOptions = { "--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=10gb --cpus-per-task=4" }
    module 'star'

    script:
    """
    STAR --runThreadN 4 \
--readFilesCommand zcat \
--genomeDir ${baseDir}/data/references/ \
--genomeFastaFiles ${genome_fna} \
--readFilesIn ${r1} \
--outFileNamePrefix ${r1.simpleName}. \
--outSAMattributes Standard \
--outSAMtype BAM SortedByCoordinate
    """

}

/*
 * Assembles the transcript using "stringtie"
 * and publish the transcript output files into the `results` folder
 */
process makeTranscript {

    input:
      file(bam_file) from bam_files
    
    output:
      file '*' into transcripts
  
    clusterOptions = { "--account=icbrbi --qos=icbrbi --time=4:00:00 --mem-per-cpu=3gb --cpus-per-task=2 --no-requeue" }
    module 'stringtie'

    """
    stringtie -e -B -p ${task.cpus} -j 2 -o ${bam_file.simpleName}_transcripts.gtf -G ${genome_gtf} ${bam_file}
    """
}

workflow.onComplete {
    def subject = 'My pipeline execution'
    def recipient = 'lboat@ufl.edu'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

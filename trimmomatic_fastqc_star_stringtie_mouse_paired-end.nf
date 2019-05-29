#!/usr/bin/env nextflow
 
/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/*_{1,2}.fq.gz"
params.genome = "/ufrc/data/reference/star/mm10_149"
params.gff = "/ufrc/data/reference/star/mm10.gtf"
   
/*
 * The reference genome file
 */
genome_file = file(params.genome)
genome_gff = file(params.gff)

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

/*
 * Step 1. QC reads
 */
process preFastqc {
    input:
    set dataset, file(reads) from fastqc_read_pairs

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    executor 'slurm'
    cpus 2
    clusterOptions = { "--account=icbrbi --qos=icbrbi-b --time=4:00:00 --mem-per-cpu=2gb --cpus-per-task=2" }
    module 'fastqc'

    script:
    """
    fastqc -q $reads
    """
}


/*
 * Step 2. Clean the reads
 */
process cleanReads {
    input:
      set dataset_id, file(forward), file(reverse) from read_pairs
 
    output:
      set dataset_id, file("${dataset_id}_R1.TRIMMED.PAIRED.fq.gz"), file("${dataset_id}_R2.TRIMMED.PAIRED.fq.gz") into paired_fastq
      set dataset_id, file("${dataset_id}_R1.TRIMMED.PAIRED.fq.gz"), file("${dataset_id}_R2.TRIMMED.PAIRED.fq.gz") into fastqc_paired_fastq
      set dataset_id, file("${dataset_id}_R1.TRIMMED.UNPAIRED.fq.gz"), file("${dataset_id}_R2.TRIMMED.UNPAIRED.fq.gz") into unpaired_fastq
 
    executor 'slurm'
    cpus 4
    clusterOptions = { "--account=icbrbi --qos=icbrbi-b --time=4:00:00 --mem-per-cpu=5gb --cpus-per-task=4" }
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
 * Step 3. QC trimmed reads
 */
process postFastqc {
    input:
    set dataset, file(forward), file(reverse) from fastqc_paired_fastq

    output:
    file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    executor 'slurm'
    cpus 2
    clusterOptions = { "--account=icbrbi --qos=icbrbi-b --time=4:00:00 --mem-per-cpu=2gb --cpus-per-task=2" }
    module 'fastqc'

    script:
    """
    fastqc -q ${forward} ${reverse}
    """
}

 
/*
 * Step 4. Maps each read-pair by using STAR
 */
process mapping {    
    input:
      file genome from genome_file
      set dataset_id, file(forward), file(reverse) from paired_fastq
  
    output:
      set dataset_id, "${dataset_id}_Aligned.sortedByCoord.out.bam" into bam_files
  
    executor 'slurm'
    cpus 4
    clusterOptions = { "--account=icbrbi --qos=icbrbi-b --time=4:00:00 --mem-per-cpu=10gb --cpus-per-task=4" }
    module 'star'

    """
    STAR \
  --genomeDir ${genome} \
  --readFilesIn ${forward},${reverse} \
  --readFilesCommand 'zcat -fc' \
  --runThreadN ${task.cpus} \
  --outSAMtype BAM SortedByCoordinate \
  --outReadsUnmapped Fastx \
  --outFilterType BySJout \
  --outFilterMultimapNmax 20 \
  --outFileNamePrefix ${dataset_id}_
    """
}
  
/*
 * Step 3. Assembles the transcript using "stringtie"
 * and publish the transcript output files into the `results` folder
 */
process makeTranscript {
    publishDir "results"
     
    input:
      set dataset_id, bam_file from bam_files
      
    executor 'slurm'
    cpus 4
    clusterOptions = { "--account=icbrbi --qos=icbrbi-b --time=4:00:00 --mem-per-cpu=5gb --cpus-per-task=4 --no-requeue" }
    module 'stringtie'

    """
    stringtie -e -B -p ${task.cpus} -j 2 -o ${dataset_id}_transcripts.gtf -G ${genome_gff} ${bam_file}
    """
}

#!/usr/bin/env nextflow

/*
 * Parameters for execution
 */
params.subject = ""
params.subject_type = "nucleotide"
params.query = ""
params.blast_type = ""
   
/*
 * Read in the files
 */
subject_fasta = file(params.subject)  
query_fasta = file(params.query)

/*
 * Index subject file
 */
process indexReference {
    publishDir "$baseDir/data/references/"

    input:
        file subject_fasta

    output:
        file "${params.subject}.*" into indices

    clusterOptions = {"-l select=1:ncpus=1:mem=5gb,walltime=04:00:00"}
    module 'singularity'

    script:
    """
    singularity run -B /scratch2 ~/singularity_containers/ncbi_blast+.simg makeblastdb -in ${subject_fasta} -dbtype ${params.subject_type} -parse_seqids
    """
}

/*
 * BLAST the sequences against the database
 */
process blastSequences {
    publishDir "$baseDir/results/"

    input:
        file "${params.subject}.*" from indices
        file query_fasta

    output:
        file "blast_result"

    clusterOptions = {"-l select=1:ncpus=2:mem=3gb,walltime=05:00:00"}
    module 'singularity'
  
    script: 
    // if-else so more different parameters may be added accordingly (at a later date)
    if( params.blast_type == "blastn"){
        """
        singularity run -B /scratch2 ~/singularity_containers/ncbi_blast+.simg blastn -db $baseDir/data/references/${params.subject} -query ${query_fasta} -evalue 1e-5 -num_threads 2 -outfmt 6 > blast_result
        """}
    else if( params.blast_type == "tblastn"){
        """
        singularity run -B /scratch2 ~/singularity_containers/ncbi_blast+.simg tblastn -db $baseDir/data/references/${params.subject} -query ${query_fasta} -evalue 1e-5 -num_threads 2 -outfmt 6 > blast_result
        """}
    else if( params.blast_type == "blastp"){
        """
        singularity run -B /scratch2 ~/singularity_containers/ncbi_blast+.simg blastp -db $baseDir/data/references/${params.subject} -query ${query_fasta} -evalue 1e-5 -num_threads 2 -outfmt 6 > blast_result
        """}
    else if( params.blast_type == "blastx"){
        """
        singularity run -B /scratch2 ~/singularity_containers/ncbi_blast+.simg blastx -db $baseDir/data/references/${params.subject} -query ${query_fasta} -evalue 1e-5 -num_threads 2 -outfmt 6 > blast_result
        """}
    else if( params.blast_type == "tblastx"){
        """
        singularity run -B /scratch2 ~/singularity_containers/ncbi_blast+.simg tblastx -db $baseDir/data/references/${params.subject} -query ${query_fasta} -evalue 1e-5 -num_threads 2 -outfmt 6 > blast_result
        """}
}

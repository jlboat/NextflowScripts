#!/bin/bash

module add nextflow

nextflow run -c config.conf ncbi-blast+.nf --subject=Sbicolor_454_v3.1.1.transcript.fa --subject_type=nucl --query=Sbicolor_454_v3.1.1.protein.fa --blast_type=tblastn 

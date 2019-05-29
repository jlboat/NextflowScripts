#!/bin/bash

module load nextflow

if [ -e ./work/ ]
then
    rm -r ./work/
    rm report.html*
    rm dag.dot*
fi

nextflow run -c config.conf rnaseq_single-end.nf -with-report -with-dag flowchart.png &> nextflow.log

# nextflow run -resume -c config.conf rnaseq_single-end.nf -with-report -with-dag flowchart.png &> nextflow.log

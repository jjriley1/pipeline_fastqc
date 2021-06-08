#!/bin/bash

if [ -d fastq.dir ]; then 
    mkdir fastq.dir/good_quality_samples
    mkdir fastq.dir/poor_quality_samples

    for filename in fastqc.dir/good_quality_samples/*.fastqc; do
        basename=$(basename "$filename")
        mv fastq.dir/${basename%%.*}* fastq.dir/good_quality_samples
    done

    for filename in fastqc.dir/poor_quality_samples/*.fastqc; do
        basename=$(basename "$filename")
        mv fastq.dir/${basename%%.*}* fastq.dir/poor_quality_samples
    done
fi

if [ -d bam.dir ]; then
    mkdir bam.dir/good_quality_samples
    mkdir bam.dir/poor_quality_samples
    
    for filename in fastqc.dir/good_quality_samples/*.fastqc; do
        basename=$(basename "$filename")
        mv bam.dir/${basename%.*}* bam.dir/good_quality_samples
    done

    for filename in fastqc.dir/poor_quality_samples/*.fastqc; do
        basename=$(basename "$filename")
        mv bam.dir/${basename%.*}* bam.dir/poor_quality_samples
    done
fi    

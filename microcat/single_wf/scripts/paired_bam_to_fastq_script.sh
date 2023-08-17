#!/bin/bash

unmapped_bam_sorted_file="$1"
unmapped_r1_fastq="$2"
unmapped_r2_fastq="$3"
unmapped_fastq="$4"
threads="$5"


echo "Starting shell script"
num_paired=$(samtools view "$unmapped_bam_sorted_file" | head -n 1000 | samtools view -c -f 1);
echo "num_paired: $num_paired"

if [ "$num_paired" -eq 1000 ]; then
    echo "Detected pair end"
    samtools fastq --threads "$threads" -n "$unmapped_bam_sorted_file" -1 "$unmapped_r1_fastq" -2 "$unmapped_r2_fastq" ;
    touch "$unmapped_fastq"
elif [ "$num_paired" -eq 0 ]; then
    echo "Detected single end"
    samtools fastq --threads "$threads" -n "$unmapped_bam_sorted_file" > "$unmapped_fastq";
    touch "$unmapped_r1_fastq" ;
    touch "$unmapped_r2_fastq" ;
else
    echo "Error: Mixed paired and unpaired reads found!";
    exit 1
fi
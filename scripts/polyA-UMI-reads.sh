#!/bin/bash
# This script extracts reads with polyA sequences from a sorted BAM file
# caveats: bam file must be sorted
# caveats: 
# usage: ./polyA-UMI-reads.sh file.bam

if [ $# -ne 1 ]; then
    echo "Usage: $0 file.bam"
    exit 1
fi
bamfile=$1
prefix="${bamfile%.bam}"

if [ ! -f $bamfile.bai ]; then
    samtools index -@8 $bamfile
fi

# rm $savepref.sam
# samtools view -H $bamfile >> $savepref.sam
samtools view -@8 -d UB $bamfile '*' | grep -vP '\tUB:Z:(\t|$)' | awk -F'\t' '$10 ~ /AAAAAA|TTTTTT/') > $prefix.polyA.sam
samtools fastq -1 $prefix.polyA.1.fq -2 $prefix.polyA.2.fq $prefix.polyA.sam

# get the UMI from the UB tag
samtools view -@8 -d UB $bamfile '*' | grep -vP '\tUB:Z:(\t|$)' #FIXME:
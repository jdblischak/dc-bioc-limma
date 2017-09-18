#!/bin/bash

index=/project2/gilad/jdblischak/dc/genome/GRCm38
fq_dir=/project2/gilad/jdblischak/dc/fastq
bam_dir=/project2/gilad/jdblischak/dc/bam

mkdir -p $bam_dir

for sample in con{1..4} gf{1..4} exgf{1..4}
do
  echo "Submitting $sample"
  sbatch --mem=12G --nodes=1 --tasks-per-node=4 --partition=broadwl --job-name=$sample-align --output=slurm-$sample-align.out run-subread.R 4 $index $fq_dir/${sample}_1.fastq.gz $fq_dir/${sample}_2.fastq.gz $bam_dir/$sample.bam
done

#!/bin/bash

bam_dir=/project2/gilad/jdblischak/dc/bam
counts_dir=data/counts
exons=/project2/gilad/jdblischak/dc/genome/exons.saf

mkdir -p $counts_dir

for sample in con{1..4} gf{1..4} exgf{1..4}
do
  echo "Submitting $sample"
  sbatch --mem=4G --nodes=1 --tasks-per-node=4 --partition=broadwl --job-name=$sample-count --output=slurm-$sample-count.out run-featurecounts.R 4 $exons $bam_dir/$sample.bam $counts_dir/$sample.txt
done

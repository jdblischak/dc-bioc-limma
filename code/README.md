# Pipeline

1. Download GEO/SRA information:
    ```
    sbatch --mem=2G --partition=broadwl code/download-info.R /project2/gilad/jdblischak/dc/data
    ```
2. Download FASTQ files:
    ```
    sbatch --mem=8G --partition=broadwl code/download-fastq.R /project2/gilad/jdblischak/dc/fastq /project2/gilad/jdblischak/dc/data/GSE75816.txt
    ```
3. Download and index genome from [Ensembl ftp](ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/):
    ```
    sbatch --mem=12G --partition=broadwl code/download-genome.R /project2/gilad/jdblischak/dc/genome
    ```
4. Map reads:
    ```
    bash code/submit-subread.sh
    ```
5. Download exons:
    ```
    sbatch --mem=2G --partition=broadwl code/download-exons.R /project2/gilad/jdblischak/dc/genome/exons.saf
    ```
6. Count reads per gene:
    ```
    bash code/submit-featurecounts.sh
    ```

Scripts adapted from my past scripts:

  * dox - [download-burridge-2016-info.R](https://github.com/jdblischak/dox/blob/master/code/download-burridge-2016-info.R)
  * dox - [download-burridge-2016-fastq.R](https://github.com/jdblischak/dox/blob/master/code/download-burridge-2016-fastq.R)
  * [midway-subread-pipline](https://github.com/jdblischak/midway-subread-pipeline)

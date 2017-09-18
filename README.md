# dc-bioc-rna-seq

Lesson on RNA-seq analysis using the data from [GSE75816][]. The scripts in
`code/` process the raw data to gene counts, which are availablin in
`data/counts/`. The lesson is in `audition.Rpres`.

## Computing environment

To setup the same computing environment, follow these steps to create a conda
environment:

1. Download and install Miniconda ([instructions](https://conda.io/miniconda.html))

1. Create the conda environment "dc" using `environment.yaml`
    ```
    conda env create --file environment.yaml
    ```
1. Activate the conda environment
    ```
    source activate
    ```

## LICENSE

[CC-BY][]. See file `LICENSE` for full license text.

[CC-BY]: https://creativecommons.org/licenses/by/4.0
[GSE75816]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75816

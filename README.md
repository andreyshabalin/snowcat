# SNOWCAT: Streamlining Our Workflow: Code Alternative for TRACTOR

Snowcat is a Genome-wide association study (GWAS) pipeline with
correction for local ancestry as estimated with RFMIX.
It follows the model originally introduced by Elizabeth G. Atkinson in the TRACTOR pipeline
(see [manuscript](https://www.nature.com/articles/s41588-020-00766-y), 
[pipeline](https://github.com/Atkinson-Lab/Tractor)).

Snowcat offers a number of advantages

1. Fast, parallel code. 
    - Allows one to rerun the whole GWAS in a matter of hours on a single computer.
2. Returns more test statistics and QC metrics for each variant
    - T-statistic and p-value for each ancestry
    - F-test and p-value for the whole model
    - Minor allele frequency (MAF) and minor allele count (MAC) for each ancestry


## Details of SNOWCAT pipeline

### Input datasets

1. Genotypes to run GWAS on
    - file format: VCF as returned by an imputation pipeline 
(e.g. [Michigan](https://imputationserver.sph.umich.edu/), 
[TopMed](https://imputation.biodatacatalyst.nhlbi.nih.gov/))







# SNOWCAT: Streamlining Our Workflow: Code Alternative for TRACTOR

Snowcat is a Genome-wide association study (GWAS) pipeline with
correction for local ancestry as estimated with RFMIX.
It follows the model originally introduced by Elizabeth G. Atkinson in the TRACTOR pipeline
(see [manuscript](https://www.nature.com/articles/s41588-020-00766-y), 
[pipeline](https://github.com/Atkinson-Lab/Tractor)).

Snowcat offers a number of improvements over the TRACTOR pipeline.

1. Fast, parallel code. 
    - Allows one to rerun the whole GWAS in a matter of hours on a single computer.
2. Returns more test statistics and QC metrics for each variant
    - T-statistic and p-value for each ancestry
    - Minor allele frequency (MAF) and minor allele count (MAC) for each ancestry






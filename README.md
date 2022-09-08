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
2. Reference genotypes data
    - VCF files, typically 1000 Genomes Reference panel (available for
[hg19](http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/) and
[GRCh38](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/)
genome builds)
    - Ancestry labels for the reference samples ([link](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel))
    - Genetic map files for the target genome ([hg19 info](https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip))

### Pipeline steps

1. Prepare reference panel files

    a. Convert the reference panel VCF files

        * Convert to binary vcf format BCF
        * Exclude SNPs with low MAF
        * Keep only genotype information (remove dosage and posterior probability info)

    b. Prepare Ancestry labels in RFMIX-friendly format
    c. Prepare genetic map files in RFMIX-friendly format

2. Convert and filter genotypes

    * Convert to binary vcf format BCF
    * Exclude SNPs with low MAF and poor imputation quality
    * Keep only genotype information (remove dosage and posterior probability info)

3. Run RFMIX estimation of local ancestry

    a. Create MSP.TSV file with ancestry info for each genomic range.

4. Convert Imputed genotypes into a custom binary format for fast parallel GWAS.

    a. Convert BCFs from step 2 into VCF for processing in R
    b. Convert VCF into custom binary format.

5. Run GWAS using RFMIX output and genotypes in binary format.

    a. Run R script for every chromosome
    b. Combine results, QC, QQ-plot and Manhattan plot.










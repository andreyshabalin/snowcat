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

## Overview of SNOWCAT pipeline

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
3. Ancestry labels for the reference samples ([link](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel))
4. Genetic map files for the target genome ([hg19 info](https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip))

### Pipeline steps

1. Prepare reference panel files
   * Convert the reference panel VCF files
     * Convert to binary vcf format BCF
     * Exclude SNPs with low MAF
     * Keep only genotype information (remove dosage and posterior probability info)
   * Prepare Ancestry labels in RFMIX-friendly format
   * Prepare genetic map files in RFMIX-friendly format
2. Convert and filter genotypes
   * Convert to binary vcf format BCF
   * Exclude SNPs with low MAF and poor imputation quality
   * Keep only genotype information (remove dosage and posterior probability info)
3. Run RFMIX estimation of local ancestry
   * Create MSP.TSV file with ancestry info for each genomic range.
4. Convert Imputed genotypes into a custom binary format for fast parallel GWAS.
   * Convert BCFs from step 2 into VCF for processing in R
   * Convert VCF into custom binary format.
5. Run GWAS with correction for local ancestry.
   * Run R script for every chromosome
   * Combine results, QC, QQ-plot and Manhattan plot.

### Install Snowcat R package

To install Snowcat R package from GitHub,
run the following code within R

```r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("andreyshabalin/snowcat")
```

# Running SNOWCAT pipeline

## Input data

1. Imputed genotypes

Imputation servers typically return imputed genotypes in files named
`chr1.dose.vcf.gz` to `chr22.dose.vcf.gz`.
The code below assumes the imputed genotypes are stored in `data_vcf_by_chr` directory:

**Genotypes:** `data_vcf_by_chr/chr"$i".dose.vcf.gz`

2. 1000 Genomes reference panel genotypes

The 1000 Genomes project genotypes are freely available online ([hg19](http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/) and
[GRCh38](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/))
in VCF format.
The code below assumes the reference panel file names to be those for hg19 genome build.

**Reference:** `ref_vcf/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`

3. Ancestry information for the reference samples

The sample map file for RFMIX. 
The RFMIX [reference manual](https://github.com/slowkoni/rfmix/blob/master/MANUAL.md)
defines its format as follows: 
It is tab delimited text with two columns. The first column gives the sample name or identifier, which must match the one used in the reference VCF/BCF. The second column is a string naming a subpopulation and may contain spaces (e.g., "European", or "East African").

The code below assumes the sample map to be located in the `map` directory:

**Sample map:** `map/1KG_map.txt`

4. Genetic map file for the used genome build.

Genetic map file for RFMIX.
The RFMIX [reference manual](https://github.com/slowkoni/rfmix/blob/master/MANUAL.md)
defines its format as follows:
The 3 columns are chromosome, physical position in bp, genetic position in cM.

Historically, the source data for creation of this file has been available at
[link](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html/files/genetic_map_b37.tar.gz) which seems to be currently broken.
You may have better luck [here](https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip) or googling for it.

The code below assumes the genetic map to be located in the `map` directory:

**Genetic map:** `map/map_file.txt`

## 1. Convert the reference panel files

We convert the reference panel VCF files to BCF to achieve several goals.
* Convert from VCF to binary vcf format BCF to greatly speed up processing by RFMIX.
* Exclude SNPs with low MAF. Very low minor allele veriants are not useful for the analysis.
* Keep only genotype information (remove dosage and posterior probability info). This greatly reduces the file size and speeds up processing by RFMIX.

The converted files are saved in `ref_bcf` directory.

<details open>
<summary>Simple bash loop (may be slow)</summary>

```bash
mkdir -p ref_bcf
for chr in {1..22}; do
 bcftools annotate \
  ref_vcf/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  -x 'FORMAT' \
  -i 'AF>.001 && AF<.999' \
  -O b9 \
  -o ref_bcf/ref_chr"$chr".bcf.gz
 bcftools index ref_bcf/ref_chr"$chr".bcf.gz
done
```

</details>

<details>
<summary>Parallel bash loop (needs a beefy machine)</summary>

```bash
mkdir -p ref_bcf
parallel --linebuffer "\
 bcftools annotate \
  ref_vcf/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  -x 'FORMAT' \
  -i 'AF>.001 && AF<.999' \
  -O b9 \
  -o ref_bcf/ref_chr{}.bcf.gz && \
 bcftools index ref_bcf/ref_chr{}.bcf.gz" ::: {1..22}
```

</details>

## 2. Convert and filter genotypes

Similarly, we convert the genotype VCF files to BCF.
In addition to previous filtering criteria, we exclude variants with imputation R2 below 0.5.

The converted files are saved in `data_bcf_by_chr_GT_QC` directory.

<details open>
<summary>Simple bash loop (may be slow)</summary>

```bash
mkdir -p data_bcf_by_chr_GT_QC
for chr in {1..22}; do
 bcftools annotate \
     data_vcf_by_chr/chr"$chr".dose.vcf.gz \
  -x 'FORMAT' \
  -i 'R2>.5 & MAF>.001' \
  -O b9 \
  -o data_bcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr"$chr".bcf.gz
 bcftools index data_bcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr"$chr".bcf.gz
done
```

</details>

<details>
<summary>Parallel bash loop (needs a beefy machine)</summary>

```bash
mkdir -p data_bcf_by_chr_GT_QC
parallel --linebuffer "\
 bcftools annotate \
     data_vcf_by_chr/chr{}.dose.vcf.gz \
  -x 'FORMAT' \
  -i 'R2>.5 & MAF>.001' \
  -O b9 \
  -o data_bcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr{}.bcf.gz && \
 bcftools index data_bcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr{}.bcf.gz" ::: {1..22}
```

</details>


## 3. Run RFMIX estimation of local ancestry

We can now run RFMIX for each chromosome.
We do not parallelize running RFMIX 
as it is very a memory intensive program.

The RFMIX output is saved in `rfmix_out` directory.

```bash
mkdir -p rfmix_out
for i in {22..1}; do
 ../rfmix/rfmix \
  --query-file=data_bcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr"$i".bcf.gz \
  --reference-file=ref_bcf/ref_chr"$i".bcf.gz \
  --sample-map=map/1KG_map.txt \
  --genetic-map=map/map_file.txt \
  --chromosome="$i" \
  --n-threads=30 \
  --output-basename=rfmix_out/rfmix_chr"$i"_R2_.5_MAF_.001
done
```

**Note:** Replace 30 in `--n-threads=30` to the number of cores on the machine being used.


## 4. Convert Imputed genotypes into a custom binary format for fast parallel GWAS.

First, we convert BCF files from step 2 to VCF for processing in R.

The converted files are saved in `data_vcf_by_chr_GT_QC`.

<details open>
<summary>Simple bash loop (may be slow)</summary>

```bash
mkdir -p data_vcf_by_chr_GT_QC
for chr in {1..22}; do
 bcftools view \
     data_bcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr"$chr".bcf.gz \
  -o data_vcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr"$chr".vcf.gz \
  -O z1
done
```

</details>

<details>
<summary>Parallel bash loop (needs a beefy machine)</summary>

```bash
mkdir -p data_vcf_by_chr_GT_QC
parallel --linebuffer "\
 bcftools view \
     data_bcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr{}.bcf.gz \
  -o data_vcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr{}.vcf.gz \
  -O z1 " ::: {1..22}
```

</details>


**Next**, we convert the genotypes to a binary format for fast random access by GWAS code.

The converted files are saved in `data_bcf_by_chr_GT_QC_fm`.

<details open>
<summary>Simple R loop (may be slow)</summary>

```R
library(snowcat)
for( chr in 1:22 ){
  vcffilename = paste0("data_vcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr",chr,".vcf.gz")
  fmnameroot = paste0("data_bcf_by_chr_GT_QC_fm/chr",chr)
  convertVCFtoFilematrix(vcffilename, fmnameroot)
}
```

</details>

<details>
<summary>Parallel bash loop (needs a beefy machine)</summary>

```bash
parallel --linebuffer "\
 Rscript -e \"\
  chr = {}; \
  library(snowcat); \
  vcffilename = paste0('data_vcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr',chr,'.vcf.gz'); \
  fmnameroot = paste0('data_bcf_by_chr_GT_QC_fm/chr',chr); \
  convertVCFtoFilematrix(vcffilename, fmnameroot);\
 \"" ::: {1..22}
```

</details>


## 5. Run GWAS with correction for local ancestry.

Running GWAS for each chromosome
can be parallelized across CPU cores of the machine.








\name{run.gwas.linear}
\alias{run.gwas.linear}
\title{Run GWAS with local ancestry}
\description{
  Run linear regression GWAS
	with correction for local ancestry
  as estimated by RFMIX.



}
\usage{
run.gwas.linear(
  chr,
  fmifilename,
  mspfilename,
  outputdir,
  cvrt,
  phenotype,
  threads = parallel::detectCores(logical = FALSE))
}
\arguments{
  \item{chr}{
    Chromosome number.
  }
  \item{fmifilename}{
    Name prefix (and possibly path) of the filematrix with genotypes
    as created by \code{\link{convertVCFtoFilematrix}}.
  }
  \item{mspfilename}{
    Name (and possibly path) of the MPS file created by RFMIX.
  }
  \item{outputdir}{
    Path for the output files.
  }
  \item{cvrt}{
    A matrix of covariates to correct for in the GWAS.
  }
  \item{phenotype}{
    A vector with the phenotype for the GWAS.
  }
  \item{threads}{
    Number of parallel threads to use for the analysis.
    Should not exceed the number of CPU cores.
  }
}
\details{
  The GWAS results are saved in the \code{outputdir} directory.
  
  The filematrix \code{gwas_chr22} contains all the GWAS results in binary numeric format.
  
  The text file \code{gwas_chr22.sumstat.txt} contains all the GWAS results in text format.
}
\value{
  The function does not return anything of value.
}
\author{
  Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\examples{
\dontrun{
# Number of principal components to correct for
npcs = 20;

# Chromosome number
chr = 22;

# File name for the filematrix with genotypes
fmifilename = paste0("data_bcf_by_chr_GT_QC_fm/chr",chr);

# Estimated local ancestry from RFMIX
mspfilename = paste0("rfmix_out/rfmix_chr",chr,"_R2_.5_MAF_.001.msp.tsv")

# Output directory
outputdir = paste0("gwas/pc",npcs, "_linear");

# Data frame with covariates
pca = fread("data_plink/pca.eigenvec",
            header = TRUE,
            data.table = FALSE);
cvrt = as.matrix(pca[ paste0("PC", seq_len(npcs)) ]);
rownames(pca) = pca$IID;
rm(npcs, pca);

# Outcome variable
samples = fread(file = paste0(fmifilename,".samples"), sep = "", quote = FALSE, header = FALSE)[[1]];
tmp = substr(samples, 1, 1);
phenotype = as.numeric(tmp == tail(tmp,1))
rm(samples, tmp);

run.gwas.linear(
    chr,
    fmifilename,
    mspfilename,
    outputdir,
    cvrt,
    phenotype)
}
}

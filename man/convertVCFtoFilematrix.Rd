\name{convertVCFtoFilematrix}
\alias{convertVCFtoFilematrix}
\title{Convert VCF File with Binary Format}
\description{
    Convert VCF file with genotype calls (GT info)
    for parallel GWAS with \code{\link{...}}.
}
\usage{
convertVCFtoFilematrix(vcffilename, fmnameroot)
}
\arguments{
    \item{vcffilename}{
        Name (and possibly path) of the VCF file (gz compressed, usually).
    }
    \item{fmnameroot}{
        Name prefix (and possibly path) of the binary output file.
    }
}
\details{
    The binary format is very economical, using only one bit per sample/variant/strand.
}
\value{
    If successful, the function returns the number of variants processed.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\examples{
\dontrun{
    chr = 22L
    vcffilename = paste0("data_vcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr",chr,".vcf.gz")
    fmnameroot = paste0("data_bcf_by_chr_GT_QC_fm/chr",chr)
    
    convertVCFtoFilematrix(vcffilename, fmnameroot)
}
}
\keyword{ bam }
\keyword{ Rbam }

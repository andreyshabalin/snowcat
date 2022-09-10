convertVCFtoFilematrix = function(vcffilename, fmnameroot){
  
  # Sanity checks
  stopifnot(is.character(vcffilename));
  stopifnot(is.character(fmnameroot));
  
  stopifnot(!is.na(vcffilename));
  stopifnot(!is.na(fmnameroot));
  
  stopifnot(length(vcffilename) == 1)
  stopifnot(length(fmnameroot) == 1)
  
  if(FALSE){
    setwd("D:/tractor/hg19/");
    chr = 22L;
    vcffilename = paste0("data_vcf_by_chr_GT_QC/GT_R2_.5_MAF_.001_chr",chr,".vcf.gz");
    fmnameroot = paste0("data_bcf_by_chr_GT_QC_fm/chr",chr);
  }

  if(!file.exists(vcffilename))
    stop("File not found: ", vcffilename);
 
  dir.create(dirname(fmnameroot), showWarnings = FALSE);
  
  if( grepl("\\.gz$",vcffilename) ){
    fid = gzfile(vcffilename, "r");
  } else {
    fid =   file(vcffilename, "r");
  }
  
  
  # Skip header, read line with sample names
  repeat{
    line = readLines(fid, 1);
    if(!grepl("^##", line))
      break;
  } # line - last non-header line, contains sample names
  
  # nchar(line)
  
  # message(line)
  spt = strsplit(line, split = "\t", fixed = TRUE)[[1]];
  snames = spt[-(1:9)];
  colnms = spt[  1:9 ];
  rm(spt, line);
  
  # Save sample names
  # library(data.table)
  # fwrite(
  #   file = paste0(fmnameroot, ".samples"),
  #   x = data.frame(samples = snames),
  #   col.names = FALSE,
  #   eol = "\n");
  writeLines(
    con = paste0(fmnameroot, ".samples"),
    text = snames);
  
  # Number of samples
  N = length(snames);
  
  # Number of bytes to store a SNP
  Nraw = (N*2 + 7) %/% 8; # length(packBits(raw(2*N), "raw"))
  
  # Read block size
  step = 1024;
  
  # Counter
  Nfilled = 0;
  
  # Indices of the values
  ind = seq_len(N*2)*2L - 1L;
  
  # vector to store genotypes to compact (must be 8x bytes)
  val2pack = raw(Nraw*8);

  # Matrix of compressed genotypes to flush to filematrix
  genoslice = matrix(raw(0),0,0);
  
  # Main loop, by 1024 (step) variants
  repeat{
    # Read a block of lines
    lines = readLines(con = fid, n = 1024);

    # If nothing is read, we are done    
    if(length(lines) == 0)
      break;

    # Exclude first 9 columns
    slice0 = gsub(
      pattern = "^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t",
      replacement =  "", 
      x = lines);
    
    # senity check
    if( !all( nchar(slice0) == (N * 4L - 1L) ) )
      stop("Does VCF include more than just GT info?")
    if( substring(slice0[1],2,2) != "|" )
      stop("Genotypes must be phased");
    
    # Matrix of compressed genotypes to flush to filematrix
    if( ncol(genoslice) != length(lines) )
      genoslice = matrix(data = raw(1), nrow = Nraw, ncol = length(lines));
    
    flip = integer(length(lines));

    # convert each variant in the slice
    for( i in seq_along(lines) ){ # i = 1
      # vector with 1 for alternative allele, 0 for reference
      values = charToRaw(slice0[i])[ind] & as.raw(1);
      # length(values);
      
      # Flip SNPs to ALT allele being minor
      if(mean(values == as.raw(0))<0.5){
        flip[i] = 1;
        values = xor(values, as.raw(1));
      }
      
      # Compact the data (1 bit per allele x SNP x sample)
      val2pack[1:(N*2L)] = values;
      packed = packBits(val2pack, "raw");
      # length(packed)
      genoslice[,i] = packed;
    }
    
    if(Nfilled == 0){
      # Create an output file
      library(filematrix);
      fm = fm.create(
        filenamebase = fmnameroot,
        nrow = Nraw,
        ncol = 0,
        type = "raw",
        size = 1);
      
      finfo = file(description = paste0(fmnameroot, ".info"), open = "wb");
      cat(file = finfo, sep = "\n", paste0(c(colnms, "FLIP"), collapse = "\t"));
    }
    
    # Save genotypes from the slice to the filematrix
    fm$appendColumns(genoslice);
    
    left = substring(text = lines, first = 1, last = nchar(lines) - N * 4L);
    cat(file = finfo, sep = "\n", paste0(left, "\t", flip));
    
    Nfilled = Nfilled + length(lines);
    message("Processed ", Nfilled, " variants in ", basename(vcffilename));
  }
  close(fm)
  close(finfo);
  message("Done processing ", basename(vcffilename));
  
  return(Nfilled);
}

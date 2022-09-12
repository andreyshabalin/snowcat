run.gwas.linear = function(
    chr,
    fmifilename,
    mspfilename,
    outputdir,
    cvrt,
    phenotype,
    threads = parallel::detectCores(logical = FALSE)){
  
  message("Starting SNOWCAT linear GWAS, chr ", chr);

  # Sanity checks
  {
    message("Checking parameters");
    
    stopifnot( length(chr) == 1 );
    stopifnot( is.numeric(chr) );
  
    stopifnot( length(fmifilename) == 1 );
    stopifnot( is.character(fmifilename) );
  
    stopifnot( length(mspfilename) == 1 );
    stopifnot( is.character(mspfilename) );
    
    stopifnot( length(outputdir) == 1 );
    stopifnot( is.character(outputdir) );
    
    dir.create(outputdir, showWarnings = FALSE, recursive = TRUE);
    stopifnot( dir.exists(outputdir) );
  
    stopifnot( is.numeric(phenotype) );
    stopifnot( is.numeric(cvrt) );
    
    Nsam = length(phenotype);
    stopifnot( NROW(cvrt) == Nsam );
    
    if( !file.exists(paste0(fmifilename, ".bmat")) )
      stop("Genotype file not found.");
    
    if( !file.exists(mspfilename) )
      stop("Genotype file not found.");
    
    stopifnot( length(threads) == 1 );
    stopifnot( is.numeric(threads) );
    stopifnot( threads >= 1 );
  } # Nsam

  # BLAS precautions
  if( threads > 1 ){
    message("Fixing BLAS settings");
    
    Sys.setenv(MKL_NUM_THREADS = 1)
    Sys.setenv(NUMEXPR_NUM_THREADS = 1)
    Sys.setenv(OMP_NUM_THREADS = 1)
  }
  
  # Info
  {
    # Info on SNPs: alleles, location, imputation quality
    message("Loading Variant info");
    info = fread(paste0(fmifilename,".info"), sep = "\t", data.table = FALSE);
    colnames(info) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "FLIP");
    Nsnp = nrow(info);
    message("Total Variants: ", Nsnp);
    
    stopifnot(all( info$CHROM == chr ));
  } # info, Nsnp
  
  # Samples
  {
    message("Loading Sample names");
    # List of samples, infer case-control status
    samples = fread(file = paste0(fmifilename,".samples"), sep = "", quote = FALSE, header = FALSE)[[1]];
    stopifnot( Nsam == length(samples) );
    message("Total samples: ", Nsam);
  } # samples, Nsam
  
  # Check genotype filematrix
  {
    # filematrix with genotypes, encoded with packBits()
    # Get dimensions
    message("Checking genotype filematrix");
    fm = fm.open(fmifilename, readonly = TRUE);
    stopifnot( ncol(fm) == Nsnp );
    stopifnot( nrow(fm) == (Nsam*2 + 7) %/% 8 );
    # fmdim = dim(fm);
    close(fm);
    rm(fm);
  } # fmdim

  # Ancestries
  {
    line1 = readLines(mspfilename, n = 1);
    line2 = strsplit(line1, split = ": ", fixed = TRUE)[[1]][2];
    spt = strsplit(line2, "\t", fixed = TRUE)[[1]];
    anslst = gsub("=.*", "", spt);
    Nans = length(anslst);
    rm(line1, line2, spt)
    message("Detected ", Nans," ancestries: ", paste0(anslst, collapse = ", "));
  } # Nans, anslst
  
  # MSP
  {
    message("Loading ancestry estimates");
    msp0 = fread(mspfilename, sep = "\t", data.table = FALSE, skip = 1);
    mspinfo = msp0[1:6];
    message("Ancestry ranges: ", nrow(mspinfo));
    # Fix end of the last interval [spos, epos]
    mspinfo$epos[nrow(mspinfo)] = msp0$epos[nrow(mspinfo)] + 1L;
    msp = as.matrix(msp0[-(1:6)]);
    rm(msp0);
    
    stopifnot( ncol(msp) == 2 * Nsam );
    stopifnot(all( colnames(msp) == t(outer( samples, 0:1, FUN = paste, sep = "."))))
    stopifnot(all( mspinfo$spos[-1] == mspinfo$epos[-nrow(mspinfo)] ));
  } # msp, mspinfo

  message("SNPs without ancestry estimate, head: ", sum(info$POS < mspinfo$spos[1]));
  message("SNPs without ancestry estimate, tail: ", sum(info$POS > tail(mspinfo$epos,1)));

  
  # Orthonormalize covariates
  {
    message("Orthonormalizing covariates");
    q = qr(cbind(1, cvrt));
    # q = qr(pcs);
    if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * NCOL(cvrt) )
      stop("Colinear or zero covariates detected");
    cvrtqr = qr.Q(q);
    rm(q);
  } # cvrtqr

  
  #### Get ready for the main loop ####
  
  # Indices of the first and second allele
  {
    ind1 = seq_len(Nsam)*2L - 1L;
    ind2 = seq_len(Nsam)*2L;
  } # ind1, ind2

  # Residualize the phenotype
  cc1 = phenotype - cvrtqr %*% crossprod(cvrtqr, phenotype);
  # cc1 = cc1 / sqrt(crossprod(cc1))

  # Create output files
  {
    message("Creating output filematrix");
    
    fmofilename = paste0(outputdir, "/gwas_chr", chr);
    # sq = 0:(nans-1L);
    rn = c(
      "chr",                             # 1
      "pos",                             # 2
      paste0("AlleleTotal_",anslst),     # 2 +
      paste0("AltAlleleCnt_",anslst),    # 2 + 1 * Nans +
      paste0("MinorAlleleFreq_",anslst), # 2 + 2 * Nans +
      paste0("Nca_",anslst),             # 2 + 3 * Nans +
      paste0("Nco_",anslst),             # 2 + 4 * Nans +
      paste0("beta_", anslst),           # 2 + 5 * Nans +
      paste0("zscore_", anslst),         # 2 + 6 * Nans +
      paste0("pvalue_", anslst),         # 2 + 7 * Nans +
      "R2",                              # 3 + 8 * Nans
      "Ftest",                           # 4 + 8 * Nans
      "Pvalue");                         # 5 + 8 * Nans
    
    fmo = fm.create(filenamebase = fmofilename, ncol = Nsnp, nrow = length(rn));
    rownames(fmo) = rn;
    colnames(fmo) = info$ID;
    close(fmo);
    
    rm(rn); # sq
    if(file.exists(paste0(fmofilename,".log")))
      file.remove(paste0(fmofilename,".log"));
  } # fmofilename
  
    
  # Loop over blocks
  
  # Common info
  ## ind1, ind2
  ## nsnp, nsam
  ## cc1
  
  # Block info:
  ## MSP of the block (msplit)
  ## set of SNPs (set)
  ## Input (genotype) and output (sum stats) file matrices
  ## (???)
  

  # Prepare the parameters
  {
    message("Preparing parameters for parallelization");
    
    fi = findInterval( c(mspinfo$spos, tail(mspinfo$epos,1)+1L), info$POS, rightmost.closed = TRUE);
    parlist = vector("list", nrow(mspinfo));
  
    for( i in seq_along(parlist) ){ # i = 1L
      set = fi[i]:(fi[i+1L] - 1L);
      mspi = c(msp[i,], use.names = FALSE) + 1L;
      levels(mspi) = as.character(1:5);
      class(mspi) = "factor";
      msplit = split(seq_along(mspi), mspi, drop = FALSE);
      # if(any(table(mspi) == 0)){
      #   stop("any");
      # }
      parlist[[i]] = list(set = set, msplit = msplit, pos = info$POS[set]);
    }
    rm(i, set, mspi, fi, msplit);

    # Order slices by size
    lns = sapply(parlist, function(x){length(x$set)});
    parlist = parlist[order(lns, decreasing = TRUE)];
    rm(lns);
    
    # Progress report numbers
    lns = sapply(parlist, function(x){length(x$set)});
    cs = c(0,cumsum(lns));
    lbl = sprintf("%.2f%% - %.2f%%", cs[-length(cs)]/tail(cs,1)*100, cs[-1]/tail(cs,1)*100);
    for( i in seq_along(parlist) ){
      parlist[[i]]$i = i;
      parlist[[i]]$lbl = lbl[i];
    }
    rm(lns, i, cs, lbl);
  }


  doSlice = function(param){
    # param = parlist[[1]]
  
    # library(filematrix)
    
    cat(as.character(Sys.time()), ", Start: slice = ", param$i, ", progress: ", param$lbl,
        "\n", sep = "", file = paste0(fmofilename,".log"), append = TRUE);
  
    # Indices of SNPs to process
    set = param$set;
    # Ancestry split
    msplit = param$msplit;
    # chromosome positions
    pos = param$pos;
    # setMKLthreads(1);
  
    # Read genotypes
    fmi = fm.open(fmifilename, readonly = TRUE);
    fmibuffer = fmi[,set];
    close(fmi);
    rm(fmi);
  
    # Open output matrix
    fmo = fm.open(fmofilename, readonly = FALSE);
    # fmobuffer = fmo[,set];
    fmobuffer = matrix(0, nrow = nrow(fmo), ncol = length(set));
  
    # Get ancestry matrix
    mat = matrix(0L, nrow = Nsam*2L, ncol = Nans);
    for( k in seq_len(Nans) ){
      mat[msplit[[k]],k] = 1L;
    }
    # crossprod(mat)
    mans = mat[ind1,] + mat[ind2,];
    rm(mat);
    
    # Allele counts
    AlleleTotal = colSums(mans);
    Nca2 = crossprod(mans, phenotype);
    Nca = Nca2 / 2;
    Nco = (AlleleTotal - Nca2)/2;
    rm(Nca2);
  
    for( j in seq_along(set) ){ # j = 1L
  
      # message(j, " of ", length(set));
      
      # if(fmobuffer[1,j]>0)
      #   next;
  
      # Load the SNP info
      snpj = as.integer(rawToBits(fmibuffer[,j])[seq_len(Nsam*2L)]);
  
      AltAlleleCnt = double(Nans);
      # AltAlleleCnt = colSums(mat2);
      
      # Split the SNP by ancestry
      mat = matrix(0L, nrow = Nsam*2L, ncol = Nans);
      for( k in seq_len(Nans) ){ # k = 1
        vals = snpj[msplit[[k]]];
        # In the rare case of all minor allele in an ancestry, exclude it (as all major allele below)
        if( any(vals == 0L) ){
          mat[msplit[[k]],k] = vals;
          AltAlleleCnt[k] = sum(vals)
        }
      }
      # crossprod(mat)
      mat2 = mat[ind1,] + mat[ind2,];
      rm(mat);
      
      # mans - matrix of ancestry dummies
      # mat2 - matrix of genotypes by ancestry
      # data.frame(anslst, AlleleTotal, AltAlleleCnt, AF = round(AltAlleleCnt/AlleleTotal,3))

  
      # working matrix for regression
      matmat = cbind(mat2, mans);
      keepset = which(colSums(matmat)>0);
      matmat = matmat[,keepset[-length(keepset)]];
      keepset5 = keepset[keepset <= Nans];
  
      # ch = chol(XX, pivot = TRUE)
      ch = chol(crossprod(matmat), pivot = TRUE)
  
      dg = which(diag(ch) < 1e-10);
      if( length(dg) > 0 ){
        pivot = attr(ch, "pivot")[dg];
        keepset = keepset[-pivot];
        matmat = matmat[,-pivot, drop = FALSE];
  
        cat("Used pivot, chr = ", chr, ", SNPid = ", set[j], "\n",
            sep = "", file = paste0(fmofilename,".pivots.txt"), append = TRUE);
      }
  
      # md = glm(cc ~ 0 + matmat + pcs1, family = "binomial");
      # ms = summary(md);
      # # ms$coefficients
      #
      # beta  = ms$coefficients[seq_along(keepset5),1];
      # se    = ms$coefficients[seq_along(keepset5),2];
      # tstat = ms$coefficients[seq_along(keepset5),3];
      # pv5   = ms$coefficients[seq_along(keepset5),4];
      #
      # an = anova(md0, md, test = "LRT");
      # pvf = an$`Pr(>Chi)`[2];
      # fff = an$Deviance[2];
  
      # Run linear regression
      X = matmat - cvrtqr %*% crossprod(cvrtqr, matmat);
      XY = crossprod(X, cc1);
      XX = crossprod(X);
  
      XX1 = solve(XX);
      XX1XY = XX1 %*% XY; # beta
      # XX1XY = solve(XX, XY);
  
      # se = sqrt(c(sig2) * diag(XX1));
  
      TSS = crossprod(cc1); # TSS = 1;
      RSS = crossprod(XY, XX1XY);
      R2 = RSS/TSS
      # R2 = crossprod(XY, XX1XY);
  
      dfFull = Nsam - ncol(cvrtqr) - ncol(matmat);
      nVarTested = ncol(matmat);
  
      fff = R2 / (1 - R2) * dfFull / nVarTested;
      pvf = pf(fff, nVarTested, dfFull, lower.tail = FALSE)
  
  
  
      beta = XX1XY[seq_along(keepset5)];
      sig2 = (TSS - RSS)/(dfFull - nVarTested);
      se = sqrt(c(sig2) * diag(XX1) * ((dfFull - nVarTested) / dfFull))[seq_along(keepset5)];
      tstat = beta/se;
      pv5 = pt(abs(tstat), df = dfFull, lower.tail = FALSE)*2;
      #
      # as.matrix(rownames(fmo))
  
  
      out = double(nrow(fmo));
      out[1] = chr;
      out[2] = pos[j]; # info$POS[set[j]];
  
      out[2 + 1:Nans] = AlleleTotal;
      out[2 + 1*Nans + 1:Nans] = AltAlleleCnt;
      out[2 + 2*Nans + 1:Nans] = pmin(AltAlleleCnt/AlleleTotal, 1 - AltAlleleCnt/AlleleTotal);
      out[2 + 3*Nans + 1:Nans] = Nca;
      out[2 + 4*Nans + 1:Nans] = Nco;
      out[2 + 5*Nans + keepset5] = beta;
      out[2 + 6*Nans + keepset5] = tstat;
      out[2 + 7*Nans + keepset5] = pv5;
      out[3 + 8*Nans] = R2;
      out[4 + 8*Nans] = fff;
      out[5 + 8*Nans] = pvf;
      
    #   if(FALSE){
    # 
    #     mf = lm(phenotype ~ 0 + matmat + cvrtqr);
    #     mfs = summary(mf);
    #     m0 = lm(phenotype ~ 0 + cvrtqr);
    #     an = anova(m0, mf)
    #     summary(mf)
    # 
    #     z = data.frame(Rcoeff = mfs$coefficients[seq_along(keepset5),1],
    #                    Scoeff = beta);
    #     z$diff = z[[1]] - z[[2]];
    #     z
    #     
    #     
    #     z = data.frame(Rse = mfs$coefficients[seq_along(keepset5),2],
    #                    Sse = se)
    #     z$diff = z[[1]] - z[[2]];
    #     z
    #     
    #     z = data.frame(Rtstat = mfs$coefficients[seq_along(keepset5),3],
    #                    Ststat = tstat );
    #     z$diff = z[[1]] - z[[2]];
    #     z
    #     
    #     z = data.frame(Rpv = mfs$coefficients[seq_along(keepset5),4],
    #                    Spv = pv5 );
    #     z$diff = z[[1]] - z[[2]];
    #     z
    # 
    #     cbind(an$F[2], fff);
    #     cbind(an$`Pr(>F)`[2], pvf);
    #   }

      fmobuffer[,j] = out;
      
      data.frame(rownames(fmo), out);
    }
  
    fmo[,set] = fmobuffer;
    close(fmo);
  
    cat(as.character(Sys.time()), ", Done:  slice = ", param$i, ", progress: ", param$lbl, "\n",
        sep = "", file = paste0(fmofilename,".log"), append = TRUE);

    return(length(set));
  } # doSlice

  
  # doSlice(param = parlist[[1]])

  # Finally run the GWAS
  if( threads > 1 ){
    message("Progress can be viewed in: ", fmofilename,".log");
    
    # library(parallel);
    tic = proc.time();
    message("Starting ", threads, " parallel processes"); 
    cl = makeCluster(threads);
    
    message("Exporting variables into the precesses");
    
    exportnames = c("fmifilename","fmofilename",
                    "chr", "Nsam", "Nsnp", "Nans",
                    "cc1", "ind1", "ind2",
                    "phenotype", "cvrtqr");
    clusterExport(cl, exportnames, envir = environment());
    rm(exportnames);
    
    message("Running GWAS in parallel");
    z = clusterApplyLB(cl, parlist, doSlice);
  
    message("Finished running GWAS");
    stopCluster(cl);
    toc = proc.time();
    message("Analysis done in ", round((toc-tic)[3]/60,1), " minutes");
  } else {
    
    # message("Progress can be viewed in: ", fmofilename,".log");
    
    # library(parallel);
    tic = proc.time();
    message("Starting single thread analysis"); 
    
    for( i in seq_along(parlist) ){ # i = length(parlist);
      message("Processing slice ", i, ", progress: ", parlist[[i]]$lbl);
      doSlice(parlist[[i]]);
    }
    
    message("Finished running GWAS");
    toc = proc.time();
    message("Analysis done in ", round((toc-tic)[3]/60,1), " minutes");
  }
  
  message("Saving GWAS output in a text file");
  mat = fm.load(fmofilename);
  df = t(mat);
  out = data.frame(
    df[,1:2], 
    ID = info$ID, 
    ref = info$REF, 
    alt = info$ALT, 
    df[,-(1:2)]);
  out = out[which(df[,1] != 0),];
  fwrite(x = out, file = paste0(fmofilename,".sumstat.txt"), sep = "\t", eol = "\n");
  
  message("Finished SNOWCAT linear GWAS, chr ", chr);
  
}

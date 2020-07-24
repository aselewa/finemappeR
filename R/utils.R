#' @export
make_ranges <- function(seqname, start, end){
  return(GenomicRanges::GRanges(seqnames = seqname, ranges = IRanges::IRanges(start = start, end = end)))
}

#' @export
clean_sumstats <- function(sumstats, cols.to.keep){
  
  stopifnot(!is.null(sumstats))
  stopifnot(length(cols.to.keep) == 8)
  
  chr <- cols.to.keep[1]
  pos <- cols.to.keep[2]
  beta <- cols.to.keep[3]
  se <- cols.to.keep[4]
  a0 <- cols.to.keep[5]
  a1 <- cols.to.keep[6]
  rs <- cols.to.keep[7]
  pval <- cols.to.keep[8]
  
  # keep SNPs in 1kg
  #sumstats <- inner_join(sumstats, snps.to.keep, by=rs)
  # Extract relevant columns
  clean.sumstats <- sumstats[ ,c(chr, pos, beta, se, a0, a1, rs, pval)]
  colnames(clean.sumstats) <- c('chr','pos','beta','se','a0','a1','snp', 'pval')
  
  # drop XY chromosomes
  clean.sumstats <- clean.sumstats[!(clean.sumstats$chr %in% c("X","Y")), ]
  print('X,Y dropped')
  # make chromosomes integers
  clean.sumstats$chr <- as.integer(clean.sumstats$chr)
  
  # Compute Zscores
  zscore <- clean.sumstats$beta/clean.sumstats$se
  clean.sumstats['zscore'] <- zscore
  clean.sumstats <- clean.sumstats[!is.na(zscore),]
  print('zscore computed')
  
  # Keep SNPs only, no indels
  nucs <- c('A','C','T','G')
  bola1 <- (clean.sumstats$a0 %in% nucs) 
  bola2 <- (clean.sumstats$a1 %in% nucs)
  clean.sumstats <- clean.sumstats[bola1 & bola2,]
  
  print('indels dropped')
  # sort by chromosome and position
  clean.sumstats <- clean.sumstats[order(clean.sumstats$chr, clean.sumstats$pos), ]
  
  # drop duplicate SNPs
  chrpos <- paste0(clean.sumstats$chr, '_', clean.sumstats$pos)
  clean.sumstats <- clean.sumstats[!duplicated(chrpos), ]
  print('duplicate snps removed')
  
  return(clean.sumstats)
}

# Assigns each SNP to one ld-block
#' @export
assign.locus.snp <- function(cleaned.sumstats, ld){
  
  ldRanges <- make_ranges(ld$X1, ld$X2, ld$X3)
  ldRanges <- plyranges::mutate(ldRanges, locus=ld$X4)
  
  snpRanges <- make_ranges(seqname = cleaned.sumstats$chr, 
                           start = cleaned.sumstats$pos, 
                           end = cleaned.sumstats$pos)
  
  snpRanges <- plyranges::mutate(snpRanges, snp=cleaned.sumstats$snp)
  
  snp.ld.overlap <- plyranges::join_overlap_inner(snpRanges, ldRanges)
  snp.ld.block <- as_tibble(snp.ld.overlap@elementMetadata)
  snp.ld.block <- snp.ld.block[!duplicated(snp.ld.block$snp), ] # some SNPs are in multiple ld-blocks due to edge of ld blocks
  cleaned.annot.sumstats <- dplyr::inner_join(cleaned.sumstats, snp.ld.block, 'snp')
  
  return(cleaned.annot.sumstats)
}


# Each annotation gets assigned SNPs based on overlap
#' @export
annotator <- function(gwas, annotations){
  
  snpRanges <- make_ranges(gwas$chr, gwas$pos, gwas$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp=gwas$snp)
  
  for(f in annotations){
    
    name <- paste0(basename(f),'_d')
    curr <- rtracklayer::import(f, format='bed')
    subdf <- IRanges::subsetByOverlaps(snpRanges, curr)
    snpsIn <- unique(subdf$snp)
    
    gwas <- dplyr::mutate(gwas, !!name := ifelse(snp %in% snpsIn,1,0))
  }
  return(gwas)
}

# Annotations for causal SNPs (apply these after fine-mapping!)
#' @export
annotator_merged <- function(gwas, annotations){
  
  snpRanges <- make_ranges(gwas$chr, gwas$pos, gwas$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp=gwas$snp)
  gwas['annots'] <- ''
  
  for(f in annotations){
    
    curr <- rtracklayer::import(f, format='bed')
    subdf <- IRanges::subsetByOverlaps(snpRanges, curr)
    snpsIn <- unique(subdf$snp)
    
    if(length(snpsIn)>0){
      curr <- gwas %>% pull(annots)
      curr <- curr[gwas$snp %in% snpsIn]
      delims <- rep(';', length(curr))
      delims[which(curr == '')] <- ''
      gwas[gwas$snp %in% snpsIn,"annots"] <- paste0(curr,delims,gsub(pattern = '.bed',replacement = '', x = basename(f)))
    }
  }
  return(gwas)
}
#' @export
merge.bigsnp.gwas <- function(gwas, bigSNP){
  
  map <- bigSNP$map
  snp_info <- map[,c('chromosome','physical.pos','allele1','allele2')]
  colnames(snp_info) <- c('chr','pos','a0','a1')
  
  matched.gwas <- as_tibble(bigsnpr::snp_match(gwas, 
                                               snp_info, 
                                               strand_flip = T, 
                                               match.min.prop = 0.1)) %>% dplyr::rename(og_index = `_NUM_ID_.ss`) %>% dplyr::rename(bigSNP_index = `_NUM_ID_`) %>% mutate(zscore = beta/se)
  
  return(matched.gwas)
}


# SUSIE related functions
#' @export
run.susie <- function(sumstats, bigSNP, ldchunk, L, prior){
  
  sub.sumstats <- sumstats[sumstats$locus == ldchunk, ]
  if(nrow(sub.sumstats) > 1){
    X <- bigSNP$genotypes[ , sub.sumstats$bigSNP_index]
    X <- scale(X, center = T, scale = T)
    zhat <- sub.sumstats$zscore
    R <- cov2cor((crossprod(X) + tcrossprod(zhat))/nrow(X))
    if(prior){
      res <- suppressWarnings(susieR::susie_rss(z = zhat, prior_weights = sub.sumstats$torus_pip, R = R, L = L, verbose = F))
    }
    else{
      res <- suppressWarnings(susieR::susie_rss(z = zhat, R = R, L = L, verbose = F))
    }
    return(res)
  }
}

# merges susie results with original summary statistics data frame
# ASSUMES L = 1! ONLY ONE CREDIBLE SET PER LOCUS!
#' @export
merge_susie_sumstats <- function(susie_results, sumstats){
  
  sumstats$susie_pip <- 0
  sumstats$CS <- 0
  loci <- names(susie_results)
  
  for(l in loci){
    n.snps <- length(susie_results[[l]]$pip)
    sumstats[sumstats$locus == as.numeric(l), "susie_pip"] <- susie_results[[l]]$pip
    
    snps.in.cs <- rep(0, n.snps)
    if(!is.null(susie_results[[l]]$sets$cs)){
      snps.in.cs[unlist(susie_results[[l]]$sets$cs$L1)] <- 1
    }
    sumstats[sumstats$locus == as.numeric(l), "CS"] <- snps.in.cs
  }
  return(sumstats)
}




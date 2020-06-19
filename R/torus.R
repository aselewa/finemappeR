
#' @title findTorus
#' @description Finds torus if it has been added to the PATH variable
#' @param path String denoting the full path to torus if known
#' @return String; torus path
findTorus <- function(path=NULL){
  if(!is.null(path)){
    TORUS <<- path
  }
  else{
    path <- system("which torus", intern=TRUE)
    if(length(path)==0){
      stop('Torus could not be found and you did not provide a path to it. Make sure Torus is installed and it is in your PATH.')
    }
    TORUS <<- path
  }
}

#' @title PrepareTorusFiles
#' @description Prepares two files necessary for running torus: z-score file and annotations file
#' @param cleaned_sumstats cleaned summary statistics from RunCleaner or other
#' @param bed_annotation a directory containing annotations in bed format. The bed file must have three columns: chr, start, end. Chromosomes should be numeric (no "chr")
#' and should be in hg19/b37 format. You will get wrong results if you use hg38 or other (reference panel is hg19).
#' @return NULL (results are written to disk in ./.temp)
PrepareTorusFiles <- function(cleaned_sumstats, bed_annotations){
  if(!exists("TORUS")){
    stop('Please run declareGlobs() first.')
  }
  stopifnot(dir.exists(bed_annotations))
  system('mkdir -p .temp')
  
  annotations <- list.files(path = bed_annotations, pattern = '*.bed', full.names = T)
  
  if(length(annotations) == 0){
    stop('No bed files/annotations were found in this directory. Are you sure this directory contains bed files?')
  }
  
  print('Annotating Snps..')
  cleaned.gwas.annots <- annotator(cleaned_sumstats, annotations = annotations)
  
  print('Writing files to temporary location..')
  data.table::fwrite(x = cleaned.gwas.annots[,-c(1:6,8:12)], file = '.temp/torus_annotations.txt.gz', quote = F, sep = '\t', col.names = T, row.names = F)
  data.table::fwrite(x = cleaned_sumstats[,c('snp','locus','zscore')], file = '.temp/torus_zscores.txt.gz', quote = F, sep = '\t', col.names = T, row.names = F)
  
  print('Done.')
}

#' @title RunTorus
#' @description executes Torus
#' @return list of enrichment results and PIP of each SNP (both tibbles)
RunTorus <- function(){
  
  if(!dir.exists('.temp')){
    stop('Cannot find annotation files. Did you run PrepareTorusfiles?')
  }
  
  args <- c('-d',
            '.temp/torus_zscores.txt.gz', 
            '-annot', 
            '.temp/torus_annotations.txt.gz',
            '--load_zval',
            '-dump_prior',
            'prior')
  
  res <- processx::run(command = TORUS, args = args, echo_cmd = TRUE, echo = TRUE)
  enrich <- as_tibble(read.table(file = textConnection(res$stdout),skip=1,header=F,stringsAsFactors = F))
  colnames(enrich) <- c("term", "estimate", "low", "high")
  
  files <- list.files(path = 'prior/', pattern = '*.prior', full.names = T)
  res2 <- processx::run(command = 'cat', args=files)
  snp_pip <- as_tibble(read.table(file = textConnection(res2$stdout),skip=1,header=F,stringsAsFactors = F))
  colnames(snp_pip) <- c("snp","torus_pip")
  
  system('rm -rf prior/')
  
  return(list(enrich=enrich, snp_pip=snp_pip))
  
}

#' @title RunTorusFDR
#' @description Runs Torus in FDR mode to estimate the FDR of each chunk containing a causal variant
#' @return tibble containing fdr of each LD block/chunk
RunTorusFDR <- function(){
  
  if(!dir.exists('.temp')){
    stop('Cannot find annotation files. Did you run PrepareTorusfiles?')
  }
  
  args <- c('-d',
            '.temp/torus_zscores.txt.gz', 
            '-annot', 
            '.temp/torus_annotations.txt.gz',
            '--load_zval',
            '-qtl')
  
  res <- processx::run(command = TORUS, args = args, echo_cmd = TRUE, echo = TRUE)
  torus_fdr <- as_tibble(read.table(file = textConnection(res$stdout),header=F,stringsAsFactors = F))
  colnames(torus_fdr) <- c("rej","region_id","fdr","decision")

  return(torus_fdr)
}


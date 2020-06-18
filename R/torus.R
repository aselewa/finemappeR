

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

findTorus('/project2/xinhe/software/dap/torus_src/torus')

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
  enrich <- read.table(file = textConnection(res$stdout),skip=1,header=F,stringsAsFactors = F)
  colnames(enrich) <- c("term", "estimate", "low", "high")
  
  res2 <- processx::run(command = 'cat prior/*')
  snp_pip <- read.table(file = textConnection(res2$stdout),skip=1,header=F,stringsAsFactors = F)
  colnames(snp_pip) <- c("SNP","PIP")
  
  return(list(enrich=enrich, snp_pip=snp_pip))
  
}

RunTorusFDR <- function(){
  
  if(!dir.exists('.temp')){
    stop('Cannot find annotation files. Did you run PrepareTorusfiles?')
  }
  
  args <- c('-d',
            '.temp/torus_zscores.txt.gz', 
            '-annot', 
            '.temp/torus_annotations.txt.gz',
            '--load_zval',
            '-qtl',
            '.temp/qtl_file.txt')
  
  res <- processx::run(command = TORUS, args = paste0(args, collapse = ' '), echo_cmd = TRUE, echo = TRUE)
  fdr_res <- readr::read_tsv(qtl_file,col_names = c("rej","region_id","fdr","decision"))
  return(fdr_res)
}

CleanTorusOutputs <- function(){
  system('rm -rf .temp')
}

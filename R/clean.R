


RunCleaner <- function(sumstats, ColsToKeep){
  if(!exists("bigSNP1kg")){
    stop('Please run declareGlobs() first.')
  }
  print('Loading summary statistics...')
  
  sumstats <- vroom::vroom(sumstats, col_names = TRUE)
  
  print('Cleaning summary statistics..')
  
  cleaned_sumstats <- clean_sumstats(sumstats, ColsToKeep)
  
  print('Matching to reference panel...')
  cleaned_sumstats <- merge.bigsnp.gwas(cleaned_sumstats, bigSNP = bigSNP1kg)
  
  print('Assining SNPs to LD blocks...')
  data('Euro_LD_Chunks', package='finemappeR')
  cleaned_sumstats <- assign.locus.snp(cleaned.sumstats = cleaned_sumstats, ld = LD_Blocks)
  
  print('Complete.')
  
  return(cleaned_sumstats)

}


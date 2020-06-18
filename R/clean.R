
bigsnp.1kg <- '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds'
ldBlocks <- readRDS('data/Euro_LD_Chunks.df.rds') 

RunCleaner <- function(sumstats, ColsToKeep){
  
  print('Loading summary statistics...')
  
  sumstats <- vroom::vroom(sumstats, col_names = TRUE)
  
  print('Cleaning summary statistics..')
  
  cleaned_sumstats <- clean_sumstats(sumstats, ColsToKeep)
  
  print('Matching to reference panel...')
  bigsnp.1kg <- snp_attach(rdsfile = bigsnp.1kg)
  cleaned_sumstats <- merge.bigsnp.gwas(cleaned_sumstats, bigsnp.1kg)
  
  print('Assining SNPs to LD blocks...')
  cleaned_sumstats <- assign.locus.snp(cleaned.sumstats = cleaned_sumstats, ld = ldBlocks)
  
  print('Complete.')
  
  return(cleaned_sumstats)

}


#' @title RunCleaner
#' @description Cleans GWAS summary statistics and adds metadata
#' @param sumstats a tibble or data frame containing raw summary statistics. Coordinates should be hg39/b37; the pipeline does not support hg38.
#' @param ColsToKeep character vector of the following columns: chr, position, allele1, allele2, beta, se, unique id, pvalue
#' @return Cleaned summary statistics + LD block of every SNP, as well as its index in the reference panel of genotypes
#' @export
RunCleaner <- function(sumstats, ColsToKeep, bigSNP){
 
  print('Loading summary statistics...')
  
  sumstats <- vroom::vroom(sumstats, col_names = TRUE)
  
  print('Cleaning summary statistics..')
  
  cleaned_sumstats <- clean_sumstats(sumstats, ColsToKeep)
  
  print('Matching to reference panel...')
  cleaned_sumstats <- merge.bigsnp.gwas(cleaned_sumstats, bigSNP = bigSNP)
  
  print('Assining SNPs to LD blocks...')
  data('Euro_LD_Chunks', package='finemappeR')
  cleaned_sumstats <- assign.locus.snp(cleaned.sumstats = cleaned_sumstats, ld = LD_Blocks)
  
  print('Complete.')
  
  return(cleaned_sumstats)

}


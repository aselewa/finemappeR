
#' @title PrepareSusieData
#' @description Adds torus results to cleaned summary statistics
#' @param sumstats a tibble or data frame containing raw summary statistics
#' @param torus_pip a tibble containing PIP of each SNP (result from RunTorus)
#' @param torus_fdr a tibble containing the FDR of each region (result from RunTorusFDR)
#' @return tibble of summary statistics updated with torus output
PrepareSusieData <- function(sumstats, torus_pip, torus_fdr, fdr_thresh=0.1){
  
  # keep loci at fdr_thresh FDR (10% by default)
  chunks <- torus_fdr$region_id[torus_fdr$fdr < fdr_thresh]
  sumstats <- sumstats[sumstats$locus %in% chunks, ]
  
  # Add Torus PIP
  sumstats <- inner_join(sumstats, torus_pip, by='snp')
  
  return(sumstats)
  
}

#' @title RunFinemapping
#' @description Runs SuSiE with L = 1
#' @param sumstats a tibble or data frame containing raw summary statistics; must have header!
#' @param priortype string; one of "torus" or "uniform"
#' @return list of finemapping results; one per LD block
RunFinemapping <- function(sumstats, priortype = "torus"){
  
  if(!exists("bigSNP1kg")){
    stop("Did you run declareGlobs()?")
  }
  
  stopifnot('torus_pip' %in% colnames(sumstats))
  chunks <- unique(sumstats$locus)
  
  if(priortype == "torus"){
    usePrior <- TRUE
    
  } 
  else if(priortype == "uniform"){
    usePrior <- FALSE
  }
  else{
    stop('prior type not recognized')
  }
  susie_res <- list()
  for(z in seq_along(chunks)){
    print(paste0('Finemapping chunks.. ',z,' of ', length(chunks)))
    susie.df <- sumstats[sumstats$locus == z, ]
    susie_res[[as.character(z)]] <- run.susie(susie.df, bigSNP1kg, z, L = 1, prior = usePrior)
  }
  
  return(susie_res)
  
}

declareGlobs <- function(bigsnp=NULL, torus=NULL){
  
  print('Declaring global variables and attaching required datasets..')
  
  if(!is.null(bigsnp)){
    bigSNP1kg <<- bigsnpr::snp_attach(bigsnp)
  } 
  else{
    bigSNP1kg <<- bigsnpr::snp_attach(rdsfile = '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds')
  }
  if(!is.null(torus)){
    TORUS <<- torus
  }
  else{
    TORUS <<- '/project2/xinhe/software/dap/torus_src/torus' 
  }
  print('Done')
}
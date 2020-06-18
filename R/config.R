
declareGlobs <- function(bigsnp=NULL, torus=NULL){
  if(!is.null(bigsnp)){
    bigsnp.1kg <<- bigsnp
  } 
  else{
    bigsnp.1kg <<- '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds'
  }
  if(!is.null(torus)){
    TORUS <<- torus
  }
  else{
    TORUS <<- '/project2/xinhe/software/dap/torus_src/torus' 
  }
}
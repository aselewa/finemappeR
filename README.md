## Introducion

finemappeR is an R package that utilizes torus and susieR to perform functionally-informed genetic finemapping. 

### Dependencies

We recommend using `R > 3.5` (has not been tested on older versions).

The package depends on TORUS, an external C++ library. To ensure this package is compiled, your system must have the following C/C++ libraries are required for compiling the source code

* GNU GSL library
* Zlib library
* Boost C++ library

### Install 

Begin by installing this repo:

```
devtools::install_github('aselewa/finemappeR')
```

After installing, check that it loads properly:

```
library(finemappeR)
```

### Reference panel

Finemapping requires linkage-disequilibrium information. We utilize the `R` package `bigsnpr` to read in PLINK files (bed/bim/fam) into `R` on the fly. If you have PLINK files, you can use `bigsnpr::snp_readBed()` to obtain a file with `.rds` extension. This RDS file is used for downstream analyses when attached with `bigsnpr::snp_attach()`.

### Vignettes

* [Quick finemapping tutorial](https://aselewa.github.io/finemappeR/articles/quick_finemapping.html)  
  
  
  

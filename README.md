## Introducion

finemappeR is an R package that utilizes torus and susieR to perform functionally-informed genetic finemapping. 

If you are at UChicago with the He Lab, installation is easy:

```
devtools::install_github('aselewa/finemappeR')
```

After installing, check that it loads properly:

```
library(finemappeR)
```

You can try to reproduce the vignette [here](vignettes/quick_finemapping.html) to ensure all functions work.

If you are an outside user, you must first install torus from here: https://github.com/xqwen/torus
Place torus on your PATH or make a note of its location for later use.

You will also need a genotype reference panel in PLINK format. We use the one readily available from 1000 Genomes Project. 
We utilize the R package `bigsnpr` to convert the PLINK files to `rds` for fast loading with `bigsnpr::snp_attach()`. 
Note the location of this `rds` file. 

With these two in hand, we are set to load the package and perform finemapping.

```
library(finemappeR)
defaultGlobs(bigSNP = path/to/bigsnp/rds, torus = path/to/torus)
```

You can try to reproduce the vignette to ensure all libraries were installed properly.

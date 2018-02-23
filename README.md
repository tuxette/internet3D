
Scripts
-------

This repository contains the scripts implementing for the network inference method of the article (Marti-Marimon *et al*, 2018). It can be installed with:

``` r
devtools::install_github("tuxette/internet3D")
```

The main functions of the package are `build.network` and `bootstrap.build`.

Supplementary files
-------------------

After installation, supplementary files with the different networks described in the article are in:

``` r
system.file("results/network0.graphml", package = "internet3D")
```

(and similarly for `network1`, `network2` and `network3`). They can be loaded within **R** with

``` r
library(igraph)
net0 <- system.file("results/network0.graphml", package = "internet3D")
net0 <- read_graph(net0, format = "graphml")
net0
```

    ## IGRAPH 99d2136 UN-- 359 2279 -- 
    ## + attr: name (v/c), id (v/c)
    ## + edges from 99d2136 (vertex names):
    ##  [1] ABCB7 --ADAMTSL3 ABCB7 --CEND1    ABCB7 --COQ7     ABCB7 --FABP3   
    ##  [5] ABCB7 --HMGB2    ABCB7 --MGST3    ABCB7 --MITF     ABCB7 --P2RX5   
    ##  [9] ABCB7 --PTMA     ABCB7 --RAB3A    ABCB7 --RBM10    ABCB7 --SELO    
    ## [13] ABI3BP--COL16A1  ABI3BP--COL5A1   ABI3BP--DPYSL3   ABI3BP--HAUS1   
    ## [17] ABI3BP--ITIH4    ABI3BP--NFATC3   ABI3BP--P4HA3    ABI3BP--RGS2    
    ## [21] ABI3BP--SDC2     ABI3BP--ZHX1     ABR   --ANXA3    ABR   --ATP1B4  
    ## [25] ABR   --CAPN10   ABR   --COL5A1   ABR   --CRLF1    ABR   --FGFR4   
    ## [29] ABR   --ISYNA1   ABR   --ITGA9    ABR   --LHX2     ABR   --NFXL1   
    ## + ... omitted several edges

The files are also directly accessible in this repository in `inst/results` and can be read with any standard tool for graph visualization (*e.g.*, Gephi, <https://gephi.org/>)

Reference
---------

Marti-Marimon M., Villa-Vialaneix N., Voillet V., Yerle-Bouissou M., Lahbib-Mansais Y., Liaubet L. (2018) A new approach of gene co-expression network inference reveals significant biological processes involved in porcine muscle development in late gestation. *Preprint*

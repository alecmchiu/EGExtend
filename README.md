# EGExtend

[![Build Status](https://travis-ci.com/alecmchiu/EGExtend.svg?branch=main)](https://travis-ci.org/alecmchiu/EGExtend)
[![codecov](https://codecov.io/gh/alecmchiu/EGExtend/branch/main/graph/badge.svg)](https://codecov.io/gh/alecmchiu/EGExtend)

An R framework for differential variance covariance testing using an eigengene-based framework. EGExtend operates by applying transformations on the data to enable standard eigengene analysis and hypothesis testing to exclusively test for variance or covariance differences as opposed to differences in means in standard eigengene analyses.

<!-- ## Citation

Please cite us at

```
Citation pending
``` -->

## Installation

EGExtend can be installed from this GitHub repository:

```r
devtools::install_github("alecmchiu/EGExtend",build_vignetttes=TRUE)
```

## Getting Started

Load EGExtend by running:

```r
library(EGExtend)
```
 
Please view the vignette for more details and examples:

 ```r
browseVignettes("EGExtend")
 ```

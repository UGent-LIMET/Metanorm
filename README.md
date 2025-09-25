# Metanorm
Robust metabolomics data normalization across scales and experimental designs

## Installing the R package

Make sure `R` (**version >= 4.4.0**) is installed on your computer:

https://cloud.r-project.org/index.html

The package can be installed directly from R as follows:

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("UGent-LIMET/Metanorm")
```

## Example usage of the R package

Load the package:

```r
library(metanorm)
```

An example dataset can be loaded from the package:
```r
load(system.file("extdata", "example.RData", package = "metanorm"))
```


The example dataset contains three objects:
- *rawdata*, a numerical matrix containing the unnormalized data
- *batch*, a vector containing for each sample run the batch to which it belongs
- *metanorm.qc*, a vector containing for each sample run whether it is a QC or another type of sample (not a QC)

Normalizing the data is achieved by calling the *metanorm* function.
```r
# normalize the first 5 compounds in the example dataset using the (default) tGAM method
normdat <- metanorm(rawdata[1:5,],       # numerical data matrix to normalize
	                model = "tGAM",      # default tGAM method
                    type = metanorm.qc,  # vector with sample types, i.e. "QC"
                                         #   and other sample types
                    QCcheck = TRUE,      # check whether QCs are representative
                    batch = batch,       # normalize by batch
                    plotdir = "~/Documents/metanormExample/")  # generate plots for
                                                               #   diagnostics
```

We can generate a PC score plot to see whether batch effects have diminished:
```r
# make a PC score plot of the data before and after normalization, label by batch
plotPCA(rawdata[1:5,], type = batch)
plotPCA(normdat[1:5,], type = batch)
```

Individual compound pre- vs. post-normalization intensity vs. order plots can be retrieved from the *plotdir* directory. These allow finegrained assessment of normalization performance and it is highly recommended to look at a decent number of these plots to judge normalization performance.

To see why this is important, you can now normalize with QC-RLSC, using, traditionally, QC samples for normalization (*QConly* argument):

```r
# normalize the first 5 compounds in the example dataset, this time using QC-RLSC
normdat2 <- metanorm(rawdata[1:5,],       # numerical data matrix to normalize
                     model = "QC-RLSC",   # use QC-RLSC
                     type = metanorm.qc,  # vector with sample types, i.e. "QC"
                                          #   and other sample types
                     QConly = TRUE,       # QC-RLSC typically used with QCs only
                     batch = batch,       # normalize by batch
                     plotdir = "~/Documents/metanormExample2/")  # generate plots for
                                                               #   diagnostics

# make a PC score plot of the data before and after normalization, label by batch
plotPCA(rawdata[1:5,], type = batch)
plotPCA(normdat2, type = batch)
```

While no apparent batch effects remain in the PC score plot, you might want to check the pre- versus post-normalization plot number 4, for example! What does the signal drift look like in the biological samples?
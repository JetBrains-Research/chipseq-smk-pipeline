if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")

BiocManager::install("Rsamtools")

if (!require(caTools))
    install.packages("caTools", repos = "http://cran.us.r-project.org")

if (!require(spp))
    install.packages("spp", dependencies=TRUE, repos = "http://cran.us.r-project.org")
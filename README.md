# SCD
  
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("sva")
    BiocManager::install("limma")
        
    source('https://raw.githubusercontent.com/jumphone/SCD/master/SCD.R')
    myscd=SCD(EXP, REF, COMBAT=TRUE)
    # Result is in "myscd$out"
    
    

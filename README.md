# SCD


    source('https://raw.githubusercontent.com/jumphone/SCD/master/SCD.R')
    
    myscd=SCD(EXP, REF, N=20, method='spearman')
    
    plot(myscd$l, col=myscd$col, pch=16)
    
    # Result is in "myscd$out"
    
    

source('SCD.R')

REF=read.table('Reference_expression.txt',sep='\t',header=T,row.names=1)
#REF[,4] = REF[,4] + REF[,5]
#REF=REF[,c(1,2,3, 4,6,7)]

REF=log(REF+1,10)
REF=apply(REF, 2, .norm_exp)
hist(REF[,1])




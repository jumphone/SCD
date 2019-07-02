


REF=readRDS('../TUMOR_REF.RDS')
.norm_exp<-function(x){
    y=x/sum(x)
    y=y*1000000
    return(y)
    }
REF=apply(REF,2,.norm_exp)

colnames(REF)=c('A','E','M','N','OL','OP','Q','NK','TA','TO','TP1','TP2','T')



###############
NUM=100
ALLR=matrix(0,ncol=NUM,nrow=13)
EXP=matrix(0,ncol=NUM,nrow=nrow(REF))

set.seed(1)
i=1
while(i<=NUM){
    this_r=runif(13)
    this_r=this_r/sum(this_r)
    this_exp=REF %*% this_r
    ALLR[,i]=this_r
    EXP[,i]=this_exp
    if(i%%10==1){print(i)}
i=i+1
}

rownames(ALLR)=colnames(REF)
rownames(EXP)=rownames(REF)
###############


#####################
getRanMap<-function(x){
  oldx=x
  x[order(x)]=sort(rnorm(length(x)))
  y=x-min(x)
  y[which(oldx==0)]=0
  return(y)
  }


set.seed(1)
REXP=apply(EXP,2,getRanMap)

REXP=REXP
  set.seed(123)
  addNOI=function(x){
     M=mean(x)
     y=x+M/3*(runif(length(x))*2-1)
     return(y)
  }

  REXP=t(apply(REXP,1,addNOI))
  REXP[which( REXP<0)]=0
  rownames( REXP)=rownames(REXP)
  colnames(REXP)=colnames(REXP)

REXP=apply(REXP,2,.norm_exp) 
################





OUT=SCD(REXP, REF, N=50, method='spearman')
plot(OUT$l)


CORMAT=cor(t(OUT$out),t(ALLR))


CORMAT=cor(OUT$exp ,EXP, method='spearman')

library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))












































.get_cor  <- function(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10, gene_check=FALSE){
    #method = "pearson", "kendall", "spearman"
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name
    
    library(parallel)
    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    
    #Step 2. calculate prob
    SINGLE <- function(i){
        library('pcaPP')
        exp_sc = as.array(exp_sc_mat[,i])
        log_p_sc_given_ref_list=rep(0,length(colname_ref))
        j=1
        while(j<=length(colname_ref)){
            exp_ref = as.array(exp_ref_mat[,j])
            log_p_sc_given_ref=cor(exp_sc,exp_ref, method=method)
            log_p_sc_given_ref_list[j]=log_p_sc_given_ref
            j=j+1}
        ################################
        #if(i%%print_step==1){print(i)}
        return(log_p_sc_given_ref_list)
        }
    #######################################
    cl= makeCluster(CPU,outfile='')
    RUN = parLapply(cl=cl,1:length(exp_sc_mat[1,]), SINGLE)
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    LOG_P_SC_GIVEN_REF = c()
    for(log_p_sc_given_ref_list in RUN){
        LOG_P_SC_GIVEN_REF=cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)}
    #######################################
    rownames(LOG_P_SC_GIVEN_REF)=colname_ref
    colnames(LOG_P_SC_GIVEN_REF)=colname_sc
    return(LOG_P_SC_GIVEN_REF)
    }



.norm_one <- function(x, cutoff=0){
    used_index=which(x>cutoff)
    non_used_index=which(x<=cutoff)
    y=x
    #y[used_index]=x[used_index]/sum(x[used_index])
    y[used_index]=scale(x[used_index]) 
    y[used_index]=pnorm(y[used_index])
    y[used_index]=y[used_index]/sum(y[used_index])
    y[non_used_index]=0
    return(y)
    }

.norm_one <- function(x){
    y=scale(x) 
    y=pnorm(y)
    y=y/sum(y)

    return(y)
    }





.scdcor  <- function(exp_sc_mat, exp_ref_mat, method='spearman'){
    #method = "pearson", "kendall", "spearman"
    method=method

    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    this_cor=cor(exp_ref_mat,exp_sc_mat,method=method)    
    return(this_cor)
    }
##






COR0=.scdcor(EXP, REF,method='spearman')
NCOR0=apply(COR0, 2, .norm_one)
cor.test(NCOR0[,1],ALLR[,1])


COR1=COR0
NCOR1=NCOR0
tmp=0

i=1
while(i<=50){
    print(i)
    EXP1=REF %*% NCOR1
    COR2=.scdcor(EXP1, REF,method='spearman')
    COR1= COR0 + COR1 - COR2
    NCOR1=apply(COR1, 2, .norm_one)
    this_exp = REF %*% NCOR1
    tmp=this_cor
    i=i+1
}



plot(NCOR1[,1],ALLR[,1])

CORMAT=cor(t(NCOR1),t(ALLR))

library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))






#####################
getRanMap<-function(x){
  oldx=x
  x[order(x)]=sort(rnorm(length(x)))
  y=x-min(x)
  y[which(oldx==0)]=0
  return(y)
  }


set.seed(1)
REXP=apply(EXP,2,getRanMap)

REXP=REXP
  set.seed(123)
  addNOI=function(x){
     M=mean(x)
     y=x+M/3*(runif(length(x))*2-1)
     return(y)
  }

  REXP=t(apply(REXP,1,addNOI))
  REXP[which( REXP<0)]=0
  rownames( REXP)=rownames(REXP)
  colnames(REXP)=colnames(REXP)

REXP=apply(REXP,2,.norm_exp) 



EXP=REXP



COR0=.scdcor(EXP, REF,method='spearman')
NCOR0=apply(COR0, 2, .norm_one)
cor.test(NCOR0[,1],ALLR[,1])


COR1=COR0
NCOR1=NCOR0

i=1
while(i<=50){
    print(i)
    EXP1=REF %*% NCOR1
    COR2=.scdcor(EXP1, REF,method='spearman')
    COR1= COR0 + COR1 - COR2
    NCOR1=apply(COR1, 2, .norm_one)
    this_exp = REF %*% NCOR1
    i=i+1
}

OUT=NCOR1
colnames(OUT)=colnames(EXP)
rownames(OUT)=colnames(REF)











COR0=.scdcor(EXP, REF,method='spearman')
NCOR0=apply(COR0, 2, .norm_one)
cor.test(NCOR0[,1],ALLR[,1])


COR1=COR0
NCOR1=NCOR0

i=1
while(i<=50){
    print(i)
    EXP1=REF %*% NCOR1
    COR2=.scdcor(EXP1, REF,method='spearman')
    COR1= COR0 + COR1 - COR2
    NCOR1=apply(COR1, 2, .norm_one)
    this_exp = REF %*% NCOR1
    i=i+1
}

OUT=NCOR1
colnames(OUT)=colnames(EXP)
rownames(OUT)=colnames(REF)

















this_out=SCD(EXP,REF,N=50)

CORMAT=cor(t(this_out),t(ALLR))

library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))


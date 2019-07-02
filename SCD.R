.norm_exp<-function(x){
    y=x/sum(x)
    y=y*1000000
    return(y)
    }

.scdcor  <- function(exp_sc_mat, exp_ref_mat, method='spearman'){
    #method = "pearson", "kendall", "spearman"
    method=method
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

.norm_one <- function(x){
    y=scale(x) 
    y=pnorm(y)
    y=y/sum(y)

    return(y)
    }


SCD <- function(EXP, REF, N=50, method='spearman'){
    ######
    REF=REF
    EXP=EXP
    method=method
    N=N
    L=c()
    
    MAT.LIST=list()
    
    
    REF=apply(REF,2,.norm_exp)
    EXP=apply(EXP,2,.norm_exp)    
    
    
    ###########
    COR0=.scdcor(EXP, REF,method=method)
    NCOR0=apply(COR0, 2, .norm_one)
    
    i=1
    print(i)
    COR1=COR0
    NCOR1=NCOR0
    colnames(NCOR1)=colnames(EXP)
    rownames(NCOR1)=colnames(REF)
    MAT.LIST=list(MAT.LIST, list(NCOR1))
    
    this_exp = REF %*% NCOR1
    cor.mat=.scdcor(this_exp ,EXP, method=method)
    L=c(L, mean(cor.mat))
    
    
    ###############
    i=2
    while(i<=N){
        print(i)
        EXP1=REF %*% NCOR1
        COR2=.scdcor(EXP1, REF, method=method)
        COR1= COR1 + COR0 - COR2
        NCOR1=apply(COR1, 2, .norm_one)
        colnames(NCOR1)=colnames(EXP)
        rownames(NCOR1)=colnames(REF)
        MAT.LIST=list(MAT.LIST, list(NCOR1))
        
        this_exp = REF %*% NCOR1
        cor.mat=.scdcor(this_exp ,EXP, method=method)
        L=c(L, mean(cor.mat))
        i=i+1
    }
    ###############
    
    max.index=which(L==max(L))
    
    
    OUT=MAT.LIST[[max.index]]
    
    EXP.OUT=REF %*% OUT
    
    RESULT=list()
    RESULT$out=OUT
    RESULT$exp=EXP.OUT
    RESULT$l=L
    RESULT$col=rep('black',length(L))
    RESULT$col[max.index]='red'
    RESULT$max.index=max.index
    return(RESULT)
    }




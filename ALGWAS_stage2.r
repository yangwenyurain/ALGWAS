# The script was used to obtain the result of ALGWAS.
# Notes:
# geno: genotype 
# phe: phenotype
# infall: the information of the genome
# infq: the information of the SNPs seleted by Adaptive Lasso
# trid: the serial number of traits
# Geno: the SNPs to be scanned
# resGeno: the SNPs that fall outside the scan window
# y1: phenotypic values to be studied
 
Pval<-function(Geno,resGeno,y1){
if(ncol(as.matrix(Geno))>1){
Genolist<-lapply(seq_len(ncol(Geno)), function(i) Geno[,i])
fit1<-lm(y1~resGeno)
rss1=sum(residuals(fit1)^2)
model<-lapply(Genolist,function(x){
fullGeno<-matrix(c(unlist(x),as.vector(resGeno)),nrow=nrow(Geno))
df1=ncol(fullGeno)-ncol(resGeno)
df2=length(y1)-ncol(fullGeno)
fit2<-lm(y1~fullGeno)
rss2=sum(residuals(fit2)^2)
Fval=(rss1-rss2)*df2/df1/rss2
pval = pf(Fval,df1,df2,lower.tail = FALSE)
return(pval)
})
} else if(ncol(as.matrix(Geno))==1){
fit1<-lm(y1~resGeno)
rss1=sum(residuals(fit1)^2)
fullGeno<-cbind(Geno,resGeno)
df1=ncol(fullGeno)-ncol(resGeno)
df2=length(y1)-ncol(fullGeno)
fit2<-lm(y1~fullGeno)
rss2=sum(residuals(fit2)^2)
Fval=(rss1-rss2)*df2/df1/rss2
pval = pf(Fval,df1,df2,lower.tail = FALSE)
model<-pval
}
return(model)
}

ALGWAS_S2<-function(geno,phe,infall,infq,trid){
y<-phe[,trid+1]
y1=y
y1[is.na(y1)]=mean(y1,na.rm=TRUE)
geno1<-as.matrix(geno[,2:ncol(geno)])
Index_snp<-matrix(0,nrow=length(infq$posx),ncol=length(infall$posx))
for (i in 1:length(infq$posx)){
   Index_snp[i,]<-(infall$posx>=infq$posx[i]-win&infall$posx<=infq$posx[i]+win)
}
# the scanned SNPs fall outside the QTL scan window
Spval<-rep(0,ncol(geno1))
covsnp<-rep(0,ncol(geno1))
S_q<-apply(as.matrix(Index_snp),2,sum)
Index_S_q<-which(S_q==0)
resGeno_0<-geno1[,colnames(geno1)%in%infq$snp]
Geno_0<-geno1[,Index_S_q]
pval0<-Pval(Geno_0,resGeno_0,y1)
Spval[Index_S_q]<-unlist(pval0)
covsnp[Index_S_q]<-nrow(infq)
# the scanned SNPs fall within the scan window of only one QTL
Index_S_q1<-which(S_q==1&S_q<length(infq$posx))
Geno_xM<-geno1[,Index_S_q1]
Index_snp1x<-Index_snp[,Index_S_q1]
Index_q1<-apply(Index_snp1x,2,function(x){which(x==1)})
Uq<-unique(Index_q1)
pvalx=list()
for (j in 1:length(Uq)){
Geno_x<-Geno_xM[,which(Index_q1==Uq[j])]
resGeno_x<-geno1[,colnames(geno1)%in%infq$snp[-Uq[j]]]
pvalx[[j]]<-Pval(Geno_x,resGeno_x,y1)
}
Spval[Index_S_q1]<-unlist(pvalx)
covsnp[Index_S_q1]<-nrow(infq)-1
# the scanned SNPs fall within the scanning windows of multiple QTLs
Index_S_q2<-which(S_q>1&S_q<length(infq$posx))
SIndex_S_q2<-length(Index_S_q2)
if (SIndex_S_q2>0){
Index_snp2<-Index_snp[,Index_S_q2]
Geno_xx<-geno1[,Index_S_q2]
pvalxx<-NULL
covsnpxx<-NULL
for (k in 1:length(Index_S_q2)){
resGeno_xx<-as.matrix(geno1[,colnames(geno1)%in%infq$snp[which(Index_snp2[,k]==0)]])
pvalxx[[k]]<-Pval(Geno_xx[,k],resGeno_xx,y1)
covsnpxx[k]<-ncol(resGeno_xx)
}
Spval[Index_S_q2]<-pvalxx
covsnp[Index_S_q2]<-covsnpxx
}
# the scanned SNPs fall within the scan window of all QTLs
SIndex_S_q3<-sum(S_q==length(infq$posx))
if (SIndex_S_q3>0){
Index_S_q3<-which(S_q==length(infq$posx))
Geno_xxx<-geno1[,Index_S_q3]
resGeno_xxx<-matrix(rep(1,nrow(geno1)),ncol=1)
pvalxxx<-Pval(Geno_xxx,resGeno_xxx,y1)
Spval[Index_S_q3]<-unlist(pvalxxx)
covsnp[Index_S_q3]<-1
}
res=data.frame(infall,covsnp=covsnp,pval=Spval)
return(res)
}
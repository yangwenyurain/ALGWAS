# The script was used to obtain the SNPs selected by Adaptive Lasso.
# Notes:
# geno: genotype 
# phe: phenotype
# infall: the information of the genome
# trid: the serial number of the trait to be studied

ALGWAS_S1<-function(geno,phe,infall,trid){
y<-phe[,trid+1]
y1=y
y1[is.na(y1)]=mean(y1,na.rm=TRUE)
fit <- msgps(as.matrix(geno[,2:ncol(geno)]),y1,penalty="alasso",gamma=1,lambda=0.001)
calasso<-coef(fit)
Nonzero<-calasso[which(calasso[,4]!=0),4]
c_snp<-names(Nonzero)
names_snp<-c_snp[2:length(c_snp)]
Index<-match(names_snp,infall$snp)
chr<-infall$chr[Index]
pos<-infall$pos[Index]
posx<-infall$posx[Index]
infq=data.frame(chr=chr,pos=pos,snp=names_snp,posx=posx,stringsAsFactors=FALSE)
return(infq)
}
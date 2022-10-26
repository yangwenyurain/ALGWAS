# The script was used to obtain the information about the genome. 
# Notes:
# geno: the genotype 
# win: the scan window of ALGWAS

Inf_geno<-function(geno,win){
names_SNP<-colnames(geno[,2:ncol(geno)])
snpall<-names_SNP
SNP_M<-unlist(strsplit(snpall,".s_"))
chrall_M<-SNP_M[seq(1,length(SNP_M),2)]
chrall<-as.numeric(unlist(strsplit(chrall_M,"r"))[seq(1,length(chrall_M<-SNP_M),2)+1])
posall<-as.numeric(SNP_M[seq(1,length(SNP_M),2)+1])
bp<-NULL
for (i in 1:max(chrall)){
bp[i]<-posall[max(which(chrall==i))]+2*win
}
chr<-c(1:max(chrall))
map<-data.frame(bp,chr)
mapx=c(0,cumsum(map$bp))
infall=data.frame(snp=snpall,chr=chrall,pos=posall,stringsAsFactors=FALSE)
infall$posx=infall$pos+mapx[infall$chr] 
return(infall)
}
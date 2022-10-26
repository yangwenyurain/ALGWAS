# ALGWAS_Demo.r

# The script was used to build a two-stage genome-wide association analysis (ALGWAS) model based on Adaptive Lasso.

# to run this ALGWAS algrithm, two files are required to be prepared as the format of the demo data:
# 1. Geno_Demo.txt: a demo genotypic data with 6000 SNP markers and 500 individuals from a maize population. The genotypic value needed to be transformed to numeric data (0/1/2).
# 2. Phe_Demo.txt: a demo phenotypic data with 500 individuals and 2 traits. The individual identifier needs to be matched with that of genotypic data.

# Several parameters and custom functions are demonstrated as followed:
# win: the scan window of ALGWAS
# trid: the serial number of the trait to be studied

rm(list=ls())
geno<-read.table("./data/Geno_Demo.txt",header=TRUE,sep="",stringsAsFactors=0,check.names=FALSE)
phe<-read.table("./data/Phe_Demo.txt",header=TRUE,stringsAsFactors=0,check.names=FALSE)

win=5*1e6
trid<-2
#trid<-1

library(msgps)
source("Inf_geno.r")
source("ALGWAS_stage1.r")
source("ALGWAS_stage2.r")

infall<-Inf_geno(geno,win)
infq<-ALGWAS_S1(geno,phe,infall,trid)
res<-ALGWAS_S2(geno,phe,infall,infq,trid)

outdir <- file.path('./res',names(phe)[trid+1])
dir.create(outdir, showWarnings = FALSE)
tr <- names(phe)[trid+1]
write.table(res,paste0(outdir,'/','res_ALGWAS_',tr,'.txt'),row.names=FALSE,quote=FALSE)



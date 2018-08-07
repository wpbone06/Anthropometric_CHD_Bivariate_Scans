#This is an R script that will do coloc analysis
#Need to do module load R first to get right version of R for coloc library
library(coloc)

args <- commandArgs(TRUE)

gwas_cad_file <- (args[1])
cad <- read.csv(gwas_cad_file, sep='\t', header=T)


eQTL_file <- (args[2])

reg <- read.csv(eQTL_file, sep='\t', header=T)
head(reg)
eQTL_N <- (args[3]) #sample size of tissue

out_file <- (args[4])

gene_name <-(args[5])

snpid_cad <- cad$rsid
snpid_reg <- reg$rsid
stopifnot( length(snpid_cad)> 1)
p_cad <- as.numeric(as.character(cad$pval))
p_reg <- as.numeric(as.character(reg$pval))
cad_N <- 34541+261984+88192+162544
#cad_N <- as.numeric(as.character(cad$Cases))+as.numeric(as.character(cad$Controls))
reg_N <- as.numeric(as.character(eQTL_N))
cad_S <- (34541+88192)/(34541+261984+88192+162544)
#cad_S <- as.numeric(as.character(cad$Cases))/cad_N
cad_maf <- as.numeric(as.character(cad$AF))
reg_maf <- as.numeric(as.character(reg$AF))
cad_maf <- pmin(cad_maf,1-cad_maf)
reg_maf <- pmin(reg_maf,1-reg_maf)


#coloc.abf.result <- coloc.abf(dataset1=list(snp=snpid_cad, pvalues=p_cad, N=cad_N, s=cad_S, type='cc'), dataset2=list(snp=snpid_reg, pvalues=p_reg, N=reg_N, type='quant'), MAF=cad_maf)
coloc.abf.result <- coloc.abf(dataset1=list(snp=snpid_cad, pvalues=p_cad, N=cad_N, s=cad_S, type='quant'), dataset2=list(snp=snpid_reg, pvalues=p_reg, N=reg_N, type='quant'), MAF=cad_maf)


write.table(coloc.abf.result$results, paste0(out_file,"results.txt"), sep='\t', row.names=F, quote=F)
write.table(coloc.abf.result$summary, paste0(out_file,"summary.txt"), sep='\t', row.names=F, quote=F)
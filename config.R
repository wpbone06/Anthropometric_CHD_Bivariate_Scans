trait1 = "CHD"
trait2 = "HDL"
trait1GWASStr = c("HDL","High-density lipoprotein","high density")
trait2GWASStr = c("Heart Disease","Coronary","Artery disease")
expPath="/project/voight_datasets/GWAS/21_CADMeta/CAD_META"
outPath="/project/voight_selscan/ksiewert/CardioMetaAnalysis/BivarAnnot/CADMeta_genvertest/HDL/jointGwasMc_HDL_formatted.txt"
trait1_BPcol = "BP"
trait2_BPcol = "hg19"
trait1_CHRcol = "CHR"
trait2_CHRcol = "SNP"
trait1_Pcol = "P-value"
trait2_Pcol = "P-value"
exp_dat = read_exposure_data("/project/voight_datasets/GWAS/21_CADMeta/CAD_META",sep="\t",snp_col="oldID",effect_allele_col="Allele1",other_allele_col="Allele2",eaf_col="Freq1",se_col="StdErr",pval_col=trait1_Pcol,beta_col="Effect")
out_dat = read_outcome_data("/project/voight_datasets/GWAS/14_lipids/jointGwasMc_HDL.txt",eaf_col="Freq.A1.1000G.EUR",effect_allele_col = "A1",other_allele_col="A2",pval_col=trait2_Pcol,samplesize_col="N",snp_col="rsid",sep="\t")
rm(list=ls())

###############################################
# Installing packages
###############################################
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")
BiocManager::install("methylGSA")
BiocManager::install("preprocessCore")
install.packages("devtools")
devtools::install_github("danbelsky/DunedinPACE")
library("ChAMP")
library("methylGSA")
library("preprocessCore")
library("DunedinPACE")

###############################################
# Setting path
# Directory must contain *.idat files and *.csv file with phenotype
###############################################
path <- "path-to-directory-with-files"
setwd(path)

###############################################
# Import and filtration
###############################################
myLoad <- champ.load(
  directory = path,
  arraytype = "EPIC",
  method = "minfi",
  methValue = "B",
  autoimpute = TRUE,
  filterDetP = TRUE,
  ProbeCutoff = 0.1,
  SampleCutoff = 0.1,
  detPcut = 0.01,
  filterBeads = TRUE,
  beadCutoff = 0.05,
  filterNoCG = TRUE,
  filterSNPs = TRUE,
  filterMultiHit = TRUE,
  filterXY = TRUE,
  force = TRUE
)

###############################################
# Functional normalization
###############################################
myNorm <- getBeta(preprocessFunnorm(myLoad$rgSet))
cpgs <- intersect(rownames(myLoad$beta), rownames(myNorm))
myNorm <- myNorm[cpgs, ]
rownames(myNorm) <- cpgs
colnames(myNorm) <- colnames(myLoad$beta)

###############################################
# Harmonization
###############################################
pd <- myLoad$pd
pd$Region <- as.factor(pd$Region)
pd$Slide <- as.factor(pd$Slide)
pd$Array <- as.factor(pd$Array)
myNorm <- champ.runCombat(
  beta = myNorm,
  pd = pd,
  variablename = "Region",
  batchname = c("Slide", "Array"),
  logitTrans = TRUE
)

###############################################
# Save
###############################################
write.csv(myNorm, file = "betas.csv")

###############################################
# DunedinPACE
###############################################
pace <- PACEProjector(myNorm)
myLoad$pd['DunedinPACE'] <- pace$DunedinPACE

###############################################
# PC clocks
# You need to setup path to 3 files from original repository (https://github.com/MorganLevineLab/PC-Clocks):
# 1) run_calcPCClocks.R
# 2) run_calcPCClocks_Accel.R
# 3) CalcAllPCClocks.RData (very big file but it is nesessary)
# You also need to apply changes from this issue: https://github.com/MorganLevineLab/PC-Clocks/issues/10
###############################################
path_pc_clocks <- "path-to-pc-clocks"
path_pc_clocks <- "D:/YandexDisk/Work/pydnameth/datasets/lists/cpgs/PC_clocks/"
source(paste(path_pc_clocks, "run_calcPCClocks.R", sep = ""))
source(paste(path_pc_clocks, "run_calcPCClocks_Accel.R", sep = ""))
myLoad$pd['Female'] <- 1
myLoad$pd[myLoad$pd$Sex == 'M', 'Female'] <- 0
pc_clocks <- calcPCClocks(
  path_to_PCClocks_directory = path_pc_clocks,
  datMeth = t(myNorm),
  datPheno = myLoad$pd
)
pc_clocks <- calcPCClocks_Accel(pc_clocks)
pc_ages <- list("PCHorvath1", "PCHorvath2", "PCHannum", "PCHannum", "PCPhenoAge", "PCGrimAge")
for (pc_age in pc_ages) {
  myLoad$pd[rownames(myLoad$pd), pc_age] <- pc_clocks[rownames(myLoad$pd), pc_age]
}

###############################################
# Save modified pheno
###############################################
write.csv(myLoad$pd, file = "pheno.csv")

###############################################
# Region-specific differences
###############################################
path_curr <- paste(path, 'region_specific', sep = '/')
dir.create(file.path(path, 'region_specific'), showWarnings = FALSE)
setwd(path_curr)

pheno <- myLoad$pd
pheno$Region <- as.factor(pheno$Region)

dmp <- champ.DMP(
  beta = myNorm,
  pheno = pheno$Region,
  compare.group = NULL,
  adjPVal = 1.0,
  adjust.method = "BH",
  arraytype = "EPIC"
)
write.csv(dmp$Central_to_Yakutia, file = "DMP.csv")
dmp_df <- data.frame(row.names(dmp$Central_to_Yakutia), dmp$Central_to_Yakutia)
colnames(dmp_df)[1] <- "CpG"
cpg_pval <- setNames(dmp_df$adj.P.Val, dmp_df$CpG)
gsea <- methylglm(
  cpg.pval = cpg_pval,
  array.type = "EPIC",
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "GO",
  minsize = 10,
  maxsize = 1000,
  parallel = TRUE
)
write.csv(gsea, file = "GO.csv", row.names=FALSE)

###############################################
# Sex-specific differences in Central
###############################################
path_curr <- paste(path, 'sex_specific_central', sep = '/')
dir.create(file.path(path, 'sex_specific_central'), showWarnings = FALSE)
setwd(path_curr)

pheno <- myLoad$pd
pheno <- pheno[pheno$Region == 'Central', ]
pheno$Sex <- as.factor(pheno$Sex)
betas_curr <- myNorm[, row.names(pheno)]

dmp <- champ.DMP(
  beta = betas_curr,
  pheno = pheno$Sex,
  compare.group = NULL,
  adjPVal = 1.0,
  adjust.method = "BH",
  arraytype = "EPIC"
)
write.csv(dmp$F_to_M, file = "DMP.csv")
dmp_df <- data.frame(row.names(dmp$F_to_M), dmp$F_to_M)
colnames(dmp_df)[1] <- "CpG"
cpg_pval <- setNames(dmp_df$adj.P.Val, dmp_df$CpG)
gsea <- methylglm(
  cpg.pval = cpg_pval,
  array.type = "EPIC",
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "GO",
  minsize = 10,
  maxsize = 1000,
  parallel = TRUE
)
write.csv(gsea, file = "GO.csv", row.names=FALSE)

###############################################
# Sex-specific differences in Yakutia
###############################################
path_curr <- paste(path, 'sex_specific_yakutia', sep = '/')
dir.create(file.path(path, 'sex_specific_yakutia'), showWarnings = FALSE)
setwd(path_curr)

pheno <- myLoad$pd
pheno <- pheno[pheno$Region == 'Yakutia', ]
pheno$Sex <- as.factor(pheno$Sex)
betas_curr <- myNorm[, row.names(pheno)]

dmp <- champ.DMP(
  beta = betas_curr,
  pheno = pheno$Sex,
  compare.group = NULL,
  adjPVal = 1.0,
  adjust.method = "BH",
  arraytype = "EPIC"
)
write.csv(dmp$F_to_M, file = "DMP.csv")
dmp_df <- data.frame(row.names(dmp$F_to_M), dmp$F_to_M)
colnames(dmp_df)[1] <- "CpG"
cpg_pval <- setNames(dmp_df$adj.P.Val, dmp_df$CpG)
gsea <- methylglm(
  cpg.pval = cpg_pval,
  array.type = "EPIC",
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "GO",
  minsize = 10,
  maxsize = 1000,
  parallel = TRUE
)
write.csv(gsea, file = "GO.csv", row.names=FALSE)

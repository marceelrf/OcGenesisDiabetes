##### Matrix ####
sampleInfo <- read.table("Data/GSE54154_Diabetes/design.csv", header=TRUE, sep=";",row.names = 1)
TMMfiles <- list.files(path = "Data/GSE54154_Diabetes/",pattern = "TMM")
Exp_list <- list()

for(i in 1:length(TMMfiles)){
  Exp_list[[i]] <- read.table(file = paste("Data/GSE54154_Diabetes/",TMMfiles[i],sep = ""),header = T,sep = "\t")
}

TMMnames <- gsub(pattern = ".data.txt",replacement = "",x = TMMfiles)
names(Exp_list) <- TMMnames

ert <- row.names(sampleInfo)
ert <- rev(ert)
ExpMatrix <- data.frame(Exp_list$GSM1308939_BMMCdb.TMM$symbol,
                        Exp_list$GSM1308939_BMMCdb.TMM$BMMCdb,
                        Exp_list$GSM1308940_BMMCdbdb.TMM$BMMCdbdb,
                        Exp_list$GSM1308941_BMGMdb.TMM$BMGMdb,
                        Exp_list$GSM1308942_BMGMdbdb.TMM$BMGMdbdb)

colnames(ExpMatrix) <- c("symbol",ert)



ExpMatrixRM <- ExpMatrix[!duplicated(ExpMatrix$symbol),]

row.names(ExpMatrixRM) <- ExpMatrixRM$symbol
ExpMatrixRM <- ExpMatrixRM[,-1]
ExpMatrixRM <- as.matrix(ExpMatrixRM)

#edgeR
dgeFULL <- DGEList(counts = ExpMatrixRM,
                   group=sampleInfo$condition)

dgeFULL <- calcNormFactors(dgeFULL)
dgeFULL <- estimateGLMCommonDisp(dgeFULL, method = "deviance", robust = T, subset = NULL)

group <- sampleInfo$Condition
design <- model.matrix(~0+group)
design
colnames(design) <- levels(group)
colnames(design)


my.contrasts <- makeContrasts(
  BMMC_db_vs_BMGM_db = BMMC_db - BMGM_db, #M2 vs M1 - Controle
  BMMC_2_vs_BMGM_2 = BMMC_2 - BMGM_2, #M2 vs M1 - Diabetes tipo 2
  BMGM_db_vs_BMGM_2 = BMGM_db - BMGM_2,#M1 controle vs M1 diabetes tipo 2
  BMMC_db_vs_BMMC_2 = BMMC_db - BMMC_2,#M2 controle vs M2 diabetes tipo 2
  levels=design)

fit <- glmFit(dgeFULL, design)
de1 <- glmLRT(fit, contrast = my.contrasts[, "BMMC_db_vs_BMGM_db"])
de2 <- glmLRT(fit, contrast = my.contrasts[, "BMMC_2_vs_BMGM_2"])
de3 <- glmLRT(fit, contrast = my.contrasts[, "BMGM_db_vs_BMGM_2"])
de4 <- glmLRT(fit, contrast = my.contrasts[,"BMMC_db_vs_BMMC_2"])


DEGs <- list()
DEGs[[1]] <- de1$table
DEGs[[2]] <- de2$table
DEGs[[3]] <- de3$table
DEGs[[4]] <- de4$table

names(DEGs) <- c("Controles","Celular_Diabetes","Diabetes_M1","Diabetes_M2")

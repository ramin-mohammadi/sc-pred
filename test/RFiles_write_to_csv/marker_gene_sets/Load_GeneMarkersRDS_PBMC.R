setwd("C:/RaminMohammadi/CelliD_R/Mostafa_Datasets")

# in right pane, then click on the .rds file and load it, then you just write them to .csv files

write.csv(PBMC_MG_10[["B"]], "PBMC_MG_10_B.csv")

write.csv(PBMC_MG_10[["CD4_T"]], "PBMC_MG_10_CD4_T.csv")

write.csv(PBMC_MG_10[["CD8_T"]], "PBMC_MG_10_CD8_T.csv")

write.csv(PBMC_MG_10[["CD14_Mono"]], "PBMC_MG_10_CD14_Mono.csv")

write.csv(PBMC_MG_10[["CD16_Mono"]], "PBMC_MG_10_CD16_Mono.csv")

write.csv(PBMC_MG_10[["DC"]], "PBMC_MG_10_DC.csv")

write.csv(PBMC_MG_10[["NK"]], "PBMC_MG_10_NK.csv")

write.csv(PBMC_MG_10[["pDCs"]], "PBMC_MG_10_pDCs.csv")


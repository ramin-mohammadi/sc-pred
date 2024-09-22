setwd("C:/RaminMohammadi/CelliD_R/Mostafa_Datasets")

# in right pane, then click on the .rds file and load it, then you just write them to .csv files

write.csv(Tabula_Lung_Total7[["alveolar"]]
          , "Tabula_Lung_Total7_alveolar.csv")

write.csv(Tabula_Lung_Total7[["B cell"]]
          , "Tabula_Lung_Total7_Bcell.csv")

write.csv(Tabula_Lung_Total7[["ciliated"]]
          , "Tabula_Lung_Total7_ciliated.csv")

write.csv(Tabula_Lung_Total7[["c-monocyte"]]
          , "Tabula_Lung_Total7_cmonocyte.csv")

write.csv(Tabula_Lung_Total7[["leukocyte"]]
          , "Tabula_Lung_Total7_leukocyte.csv")

write.csv(Tabula_Lung_Total7[["endothelial"]]
          , "Tabula_Lung_Total7_endothelial.csv")


write.csv(Tabula_Lung_Total7[["mast"]]
          , "Tabula_Lung_Total7_mast.csv")

write.csv(Tabula_Lung_Total7[["myeloid"]]
          , "Tabula_Lung_Total7_myeloid.csv")

write.csv(Tabula_Lung_Total7[["NK"]]
          , "Tabula_Lung_Total7_NK.csv")

write.csv(Tabula_Lung_Total7[["nc-monocyte"]]
          , "Tabula_Lung_Total7_nc_monocyte.csv")

write.csv(Tabula_Lung_Total7[["stromal"]]
          , "Tabula_Lung_Total7_stromal.csv")

write.csv(Tabula_Lung_Total7[["T cell"]]
          , "Tabula_Lung_Total7_Tcell.csv")

write.csv(Tabula_Lung_Total7[["pneumocyte-II"]]
          , "Tabula_Lung_Total7_pneumocyteII.csv")

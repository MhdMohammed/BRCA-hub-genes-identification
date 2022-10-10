

#https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/casestudy.html

library(SummarizedExperiment)
library(TCGAbiolinks)
query.exp <- GDCquery(project = "TCGA-BRCA", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))
GDCdownload(query.exp)
BRCA.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "BRCAExp.rda")

# get subtype information
dataSubt <- TCGAquery_subtype(tumor = "BRCA")

# get clinical data
dataClin <- GDCquery_clinic(project = "TCGA-BRCA","clinical") 

# Which samples are primary solid tumor
dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP") 
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")

dataPrep <- TCGAanalyze_Preprocessing(object = BRCA.exp, cor.cut = 0.6)                      

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
                            mat2 = dataFilt[,dataSmTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")



ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",
                                RegulonList = rownames(dataDEGs))  

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = rownames(dataDEGs),
                        nBar = 20)


group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))

dataSurv <- TCGAanalyze_SurvivalKM(clinical_patient = dataClin,
                                   dataGE = dataFilt,
                                   Genelist = rownames(dataDEGs),
                                   Survresult = FALSE,
                                   ThreshTop = 0.67,
                                   ThreshDown = 0.33,
                                   p.cut = 0.05, group1, group2)



require(dnet)  # to change
org.Hs.string <- dRDataLoader(RData = "org.Hs.string")

TabCoxNet <- TCGAvisualize_SurvivalCoxNET(dataClin,
                                          dataFilt, 
                                          Genelist = rownames(dataSurv),
                                          scoreConfidence = 700,
                                          org.Hs.string = org.Hs.string,
                                          titlePlot = "Case Study n.1 dnet")


MultiData$stages <- ifelse(MultiData$stages=="stage iiic", "stageiii", 
                    ifelse(MultiData$stages=="stage iiib", "stageiii", 
                    ifelse(MultiData$stages=="stage iiia", "stageiii", 
                    ifelse(MultiData$stages=="stage iii", "stageiii", 
                    ifelse(MultiData$stages=="stage iib", "stageii", 
                    ifelse(MultiData$stages=="stage iia", "stageii", 
                    ifelse(MultiData$stages=="stage ii", "stageii", 
                    ifelse(MultiData$stages=="stage ib", "stagei", 
                    ifelse(MultiData$stages=="stage ia", "stagei", 
                    ifelse(MultiData$stages=="stage i", "stagei", 
                    ifelse(MultiData$stages=="stage iv", "stageiv", 
                    ifelse(MultiData$stages=="stage x", "stagex", "NA"))))))))))))
MultiDataDrop <- MultiData[MultiData$stages != "NA",]

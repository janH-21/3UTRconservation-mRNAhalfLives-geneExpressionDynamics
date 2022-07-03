##### function definition ##### 
# function to export results of GO analysis as .csv files
exportGOtable <- function(resultData){
  resultData <- subset(resultData, select = c(description, 
                                              geneSet, 
                                              size, 
                                              expect, 
                                              overlap, 
                                              enrichmentRatio, 
                                              pValue, 
                                              FDR))
  colnames(resultData) <- c("Description", 
                            "GO term", 
                            "Number of genes", 
                            "Expected number of genes of interest",
                            "Detected number of genes of interes",
                            "Enrichment",
                            "p",
                            "FDR")
  return(resultData)
}

##### load packages #####
require(WebGestaltR)
require(tidyverse)
require(kableExtra)

##### GO analysis ##### 
## load data
hg38_pivot_unique <- read.csv("dataR/hg38_pivot_unique.csv")

## save reference gene set
write.table(hg38_pivot_unique$gene_name, 
            "GOanalysisR/reference.txt", 
            sep="\t", col.names = F, row.names = F, quote = FALSE)

## iterate several cutoffs for conservation and length for sensitivity analysis
for(cons in c(0.6, 0.7, 0.8, 0.9)){
  for(len in c(1000, 1500, 2000, 5000)){
    if(cons == 0.6 | len == 1000){
      
      # choose gene list for given cutoff
      geneList <- hg38_pivot_unique %>% dplyr::filter(threeUTR_length > len, threeUTR_conservation > cons)
      # save list obtained with given cutoff
      write.table(geneList$gene_name, 
                  paste("GOanalysisR/", len, "_", cons, ".txt", sep = ""), 
                  sep="\t", col.names = F, row.names = F, quote = FALSE)
      
      # run webgestalt GO analysis for categories biological process and molecular function
      for(process in c("geneontology_Biological_Process", "geneontology_Molecular_Function")){
        # run top10 analysis for first look and fdr based analysis for proper analysis
        for (method in c("top","fdr")) {
          
          pName = paste(format(Sys.time(), "%Y-%m-%d_%H-%M"), len, cons, process, method, sep="_")
          
          webGestaltResults <- WebGestaltR(enrichMethod = "ORA", organism = "hsapiens", enrichDatabase = process,
                                           interestGene = geneList$gene_name, interestGeneType = "ensembl_gene_id", 
                                           collapseMethod = "mean",
                                           referenceGene = hg38_pivot_unique$gene_name, referenceGeneType = "ensembl_gene_id",
                                           minNum = 5, maxNum = 2000, 
                                           sigMethod = method, fdrMethod = "BH", fdrThr = 0.01, topThr = 20, 
                                           reportNum = 20, perNum = 1000, gseaP = 1,
                                           isOutput = TRUE, outputDirectory = "GOanalysisR", 
                                           projectName = pName,
                                           hostName = "http://www.webgestalt.org/"
          )
          
        }
      }
    }
  }
}

##  set directory 
wd <- "GOanalysisR/"
GOdirs <- dir(path = wd, pattern = "Project.*fdr") 

## iterate results
for(GOdirPath in GOdirs){
  
  # choose current result file
  GOresultFile <- dir(path = paste(wd, GOdirPath, sep =""), pattern = "_results_")
  
  # create output file
  fileOut <- GOresultFile %>% 
    substr(20, nchar(GOresultFile)-4) %>%
    paste(wd, ., ".csv",sep="")
  
  # read, format and write result data
  read.csv(paste(wd, GOdirPath, "/", GOresultFile, sep = ""), header = TRUE, sep = "\t", dec = ".") %>% 
    exportGOtable %>%
    write.csv(fileOut)
  
}

## find entries duplicated by WebGestalt during  mapping
mappedIDs <- read.csv("GOanalysis/Project_2020_11_22_21_09_1000_0_6_geneontology_Molecular_Function_fdr/interestingID_mappingTable_2020_11_22_21_09_1000_0_6_geneontology_Molecular_Function_fdr.txt", header = TRUE, sep="\t")

doubleEntries <- mappedIDs[mappedIDs$userId%in% mappedIDs[duplicated(mappedIDs$userId),"userId"],] 
doubleEntries %>% kbl() %>% kable_styling(full_width = FALSE)  

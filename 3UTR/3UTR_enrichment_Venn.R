## load packages
require(tidyverse)
require(eulerr)

## load data
hg38_pivot_TF <- read.csv("dataR/hg38_pivot_TF.csv")

## get set sizes
TFs <- hg38_pivot_TF %>%
  dplyr::filter(isTF == TRUE) %>% 
  nrow # 1797 (old version based on hg19: )
all <- hg38_pivot_TF %>% 
  nrow # 1811 (old version based on hg19: )
lenCons <- hg38_pivot_TF %>% 
  dplyr::filter(threeUTR_length>1000, threeUTR_conservation > 0.6) %>% 
  nrow # 845 (old version based on hg19: )
lenConsTF <- hg38_pivot_TF %>% 
  dplyr::filter(threeUTR_length>1000, threeUTR_conservation > 0.6, isTF == TRUE) %>% 
  nrow # 226 (old version based on hg19: )

## enrichment ratio
rEnrichment <- ((lenConsTF/lenCons)/(TFs/all))

## hypergeometric test
pEnrichment <- phyper(lenConsTF-1,
                      TFs,
                      (all - TFs),
                      lenCons,
                      lower.tail = FALSE)
cat("Enrichment ratio = ", rEnrichment, "\n", "p (Hypergeometric) = ", pEnrichment, sep="")

## Venn diagram
abc <- c(A = all - TFs - lenCons + lenConsTF, 
         B = 0, 
         C = 0, 
         "A&B" = lenCons - lenConsTF, 
         "A&C" = TFs - lenConsTF, 
         "A&B&C" = lenConsTF)

g <- plot(euler(abc), 
          labels = c(paste("all genes,\n n = ", toString(all-lenCons-TFs+lenConsTF)),
                     paste("\n transcription factors,\n n = ", toString(TFs-lenConsTF)),
                     paste("\n \n genes with 3'UTR lenght > 1000bp und conservation >0.6,\n n = ",toString(lenCons-lenConsTF),
                           "\n\n TFs with 3'UTR length > 1000bp und conservation >0.6,\n n = ",toString(lenConsTF)
                     )),
          fill = c("orange", "blue", "green", "red"),
          alpha = c(0.6,0.7,0.6),
          edges = "transparent")
g

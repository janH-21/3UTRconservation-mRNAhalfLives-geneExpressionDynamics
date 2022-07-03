## load packages
require(tidyverse)
require(ggpubr)

## load data
hg38_TF <- read.csv("dataR/hg38_pivot_TF.csv")

## scatter plot len vs. cons
scatterPlot <- ggplot(hg38_TF, aes(x=threeUTR_conservation, y=threeUTR_length, color=isTF)) +
  geom_point(size = 0.5)

## cons histogram + density
conHist <- ggplot(hg38_TF, aes(x=threeUTR_conservation, color=isTF)) +
  geom_histogram(alpha=0, size=0.33) + xlim(0, 1)
conDens <- ggplot(hg38_TF, aes(x=threeUTR_conservation, color=isTF)) +
  geom_density(alpha=0, size=0.5) + xlim(0, 1)

## len histogram + density
lenHist <-  ggplot(hg38_TF, aes(x=threeUTR_length, color=isTF)) +
  geom_histogram(alpha=0 ,size=0.33) + xlim(0, max(hg38_TF$threeUTR_length)) + coord_flip()
lenDens <- ggplot(hg38_TF, aes(x=threeUTR_length, color=isTF)) +
  geom_density(alpha=0, size=0.5) + xlim(0, max(hg38_TF$threeUTR_length)) + coord_flip()

## arrange plots
gCons <- ggarrange(conHist + rremove("x.text") + rremove("y.text"), 
                   conDens + rremove("x.text") + rremove("y.text"), 
                   ncol = 1, legend = FALSE)
gLen <- ggarrange(lenHist + rremove("x.text") + rremove("y.text"), 
                  lenDens + rremove("x.text") + rremove("y.text"), 
                  nrow = 1, legend = FALSE)
gScatterAll <- ggarrange(scatterPlot + rremove("x.text") + rremove("y.text"),
                         gLen, gCons,
                         nrow = 2, ncol = 2, legend = FALSE)
gScatterAll


## NOTE: plot axes
# - length: 0, 10k, 20k 
# - length count: 0:2000:6000
# - cons: 0, 0.25, 0.5, 0.75, 1
# - cons count: 0:500:2000 
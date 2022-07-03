# 3UTRconservation-mRNAhalfLives-geneExpressionDynamics

1. Folder "3UTR" contains R scripts to process analyse and plot 3'-UTR data. Start with 3UTR_prepareData.R for formating and cleaning of data mapped from the UCSC genome browser, as well as to add the data from the other data sources. All other scripts in this folder can then be used with the data sets saved from 3UTR_prepareData.R in no particular order.
2. Folder "model" contains C++ programms to run the model of transcriptional dynamics. Ideally, don't run at once, but only some input frequencies at once, then adapt max. time and output sampling rate in order to reduce overall run time and output data size.
3. Folder "modelAnalysis" contains R scripts to analyse output data from the C++ models and plot the results.

Figures and statistics are directly reproduced by the scripts provided here.

#### load packages ####
require(tidyverse)
require(reshape2)
require(ensembldb)

####  read .gtf data ####
hg38 <- read.delim("dataR/hg38.ensGene.conservation_calculation.gtf", header = FALSE, sep = "\t", dec = ".")

#### format data #### 
## prepare empty vectors to hold data
hg38_gene_id <- vector("list", nrow(hg38))
hg38_transcript_id <- vector("list", nrow(hg38))
hg38_gene_name <- vector("list", nrow(hg38))

## iterate hg38, split column containing IDs into separate columns
for(iter in 1:nrow(hg38)){
  colContent <- toString(hg38$V9[iter])
  colContent <- (strsplit(colContent, " "))
  colContent <- colContent[[1]]
  for(jter in 1:length(colContent)){
    if(colContent[jter]=="gene_id"){
      hg38_gene_id[iter] <- strsplit(colContent[jter+1], ";")
    }
    if(colContent[jter]=="transcript_id"){
      hg38_transcript_id[iter] <- strsplit(colContent[jter+1], ";")
    }
    if(colContent[jter]=="gene_name"){
      hg38_gene_name[iter] <- strsplit(colContent[jter+1], ";")
    }
  }
}

## save intermediate data sets
df_hg38_gene_id <- hg38_gene_id %>% unlist %>% data.frame %>% 'colnames<-' (c("gene_id"))
write.csv(df_hg38_gene_id, "dataR/df_hg38_gene_id.csv", row.names = FALSE)

df_hg38_gene_name <- hg38_gene_name %>% unlist %>% data.frame %>% 'colnames<-' (c("gene_name"))
write.csv(df_hg38_gene_name, "dataR/df_hg38_gene_name.csv", row.names = FALSE)

df_hg38_transcript_id <- hg38_transcript_id %>% unlist %>% data.frame %>% 'colnames<-' (c("transcript_id"))
write.csv(df_hg38_transcript_id, "dataR/df_hg38_transcript_id.csv", row.names = FALSE)

## convert to DF
df_chromosome <- hg38$V1 %>% data.frame %>% 'colnames<-' (c("chromosome"))
df_sequence_type <- hg38$V3 %>% data.frame %>% 'colnames<-' (c("sequence_type"))
df_strand <- hg38$V7 %>% data.frame %>% 'colnames<-' (c("strand"))
df_length <- hg38$V5 - hg38$V4  %>% data.frame %>% 'colnames<-' (c("length"))
df_conservation <- hg38$V12 %>% data.frame %>% 'colnames<-' (c("conservation"))
hg38_ordered <- cbind(df_hg38_gene_name,
                      df_hg38_gene_id,
                      df_hg38_transcript_id,
                      df_chromosome,
                      df_strand,
                      df_sequence_type,
                      df_length,
                      df_conservation)

## save intermediate data set
write.csv(hg38_ordered, "dataR/hg38_ordered.csv", row.names = FALSE)

#### compile data by transcript #### 
## load data, if starting from hete
# hg38_ordered <- read.csv("dataR/hg38_ordered.csv")
hg38_pivot <- hg38_ordered %>%
  melt(id.vars = c("transcript_id", "sequence_type"), measure.vars = c("length", "conservation")) %>%
  dcast(transcript_id ~ sequence_type + variable, mean)

## replace column names to avoid errors due to column names starting with numbers
colnames(hg38_pivot) <-c("transcript_id",            
                         "threeUTR_length",
                         "threeUTR_conservation",
                         "fiveUTR_length",           
                         "fiveUTR_conservation",        
                         "CDS_length",               
                         "CDS_conservation",         
                         "exon_length",             
                         "exon_conservation",
                         "start_codon_length",
                         "start_codon_conservation",
                         "stop_codon_length",
                         "stop_codon_conservation",  
                         "transcript_length",
                         "transcript_conservation")

#### add gene names ####
unique(hg38_ordered$gene_id == hg38_ordered$gene_name) # gene names == gene ids? -> true
mapIDs <- dplyr::select(hg38_ordered, c(1,3)) %>% subset(!duplicated(transcript_id)) # select ENSG and ENST IDs as mapping
hg38_pivot_all <- merge(mapIDs, hg38_pivot) # actual mapping

####  select data of interest, clean NAs and sequences <10bp length #### 
hg38_pivot_cleaned <- hg38_pivot_all %>% 
  dplyr::select(c("gene_name", 
                  "transcript_id", 
                  "threeUTR_length", 
                  "threeUTR_conservation", 
                  "fiveUTR_length", 
                  "fiveUTR_conservation",        
                  "CDS_length",               
                  "CDS_conservation")) %>%
  dplyr::filter(!is.na(threeUTR_length),
                !is.na(threeUTR_conservation),
                !is.na(fiveUTR_length),
                !is.na(fiveUTR_conservation),
                !is.na(CDS_length),
                !is.na(CDS_conservation)) %>%
  dplyr::filter(threeUTR_length>10,
                fiveUTR_length>10,
                CDS_length>10)

## unique based on gene IDs
hg38_pivot_unique <- hg38_pivot_cleaned %>%
  melt(id.vars = c("gene_name"), measure.vars = c("threeUTR_length", 
                                                  "threeUTR_conservation", 
                                                  "fiveUTR_length", 
                                                  "fiveUTR_conservation",        
                                                  "CDS_length",               
                                                  "CDS_conservation")) %>%
  dcast(gene_name ~ variable, mean)

## save intermediate data set
write.csv(hg38_pivot_unique, "dataR/hg38_pivot_unique.csv", row.names = FALSE)


#### add gene symbols for WebGestalt Analysis #### 
hg38_symbl <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                keys= hg38_pivot_unique$gene_name, 
                                keytype = "GENEID", 
                                columns = c("SYMBOL","GENEID"))
hg38_pivot_symbl <- merge(hg38_pivot_unique, 
                          hg38_symbl, 
                          by.x = "gene_name", 
                          by.y = "GENEID",
                          all.x = FALSE, 
                          all.y = FALSE)

#### add TF data ####
TFgenom <- read.csv("dataR/human_tf_genomatix.txt", header = FALSE) 
colnames(TFgenom) <- "gene_name"
TFgenom$isTF <- replicate(nrow(TFgenom), TRUE)
hg38_pivot_TF <- merge(hg38_pivot_symbl, TFgenom, by.x = "SYMBOL", by.y = "gene_name", all.x = TRUE, all.y = FALSE)
hg38_pivot_TF$isTF[is.na(hg38_pivot_TF$isTF)] <- FALSE
## save intermediate data
write.csv(hg38_pivot_TF, "dataR/hg38_pivot_TF.csv", row.names = FALSE)


##### add halfLife data ##### 
halfLife <- read.csv("dataR/human_table2.txt", sep="\t")
hg38_pivot_half <- merge(hg38_pivot_TF, halfLife, by.x = "SYMBOL", by.y = "name", all.x = FALSE, all.y = FALSE)
## save intermediate data
write.csv(hg38_pivot_half, "dataR/hg38_pivot_half.csv", row.names = FALSE)

####################################
#### Script for making a dataframe containing Penstemon species and corresponding OTU taxonomy and abundances from multinomial data
#Results in two .csv files that merge OTU ID's, taxonomy, Penstemon ID, and multinomial data, one for 16S and one for ITS

#Load packages
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scico)
require(RDPutils)

#setwd
setwd("/Users/ebentley/Desktop/PenstMicrobe-master")

#### Data wrangling 16S multinomial ####
#calling both multinomial16S and tax16Sonomy and doing some formatting for convenient use later
metadat <- read.csv("./metadata_wrangled.csv", stringsAsFactors = F)
identity<- metadat[3] #use 22 when you are trying to get the location code in there.
multinomial16S <- read.csv("multinomial_estimates_16s_by_taxon.csv", stringsAsFactors = F, header = T) 
multinomial16S <-cbind(identity, multinomial16S)
multinomial16S <- data.frame(multinomial16S[,-1], row.names = multinomial16S[,1])

#standardize by the ISD
multinomial16S <- multinomial16S[, 2:length(multinomial16S)] /multinomial16S[,names(multinomial16S)== "ISD"]

#get rid of the mtDNA, Cholroplast, plant, and duds columns in the multinomial16S dataframe
multinomial16S<- multinomial16S[, 1:(length(multinomial16S)-4)]

#import the taxonomy table with 16S data 
tax16S <- import_sintax_file( "penst16S_ESV.gg.sintax", confidence = 0.8)
tax16S <- as.data.frame(tax16S)
tax16S[tax16S == "uncl_Domain"] <- "NA"

#Split the names at the _ so that we can separate out the bacterial names from the syntax16S deliniations of d,k,p,c,etc. 
#Run this twice--some OTUs have two sections of "_"s

tax16S$Domain<- gsub("^.*?_", "", tax16S$Domain)
tax16S$Phylum<- gsub("^.*?_", "", tax16S$Phylum)
tax16S$Class<- gsub("^.*?_", "", tax16S$Class)
tax16S$Order<- gsub("^.*?_", "", tax16S$Order)
tax16S$Family<- gsub("^.*?_", "", tax16S$Family)
tax16S$Genus<- gsub("^.*?_", "", tax16S$Genus)
tax16S$Species<- gsub("^.*?_", "", tax16S$Species)

#make sure no unwanted OTU's are present
tax16S<-subset(tax16S, Family != "mitochondria")
tax16S<-subset(tax16S, Class != "Chloroplast")
tax16S[tax16S == "uncl_Domain"] <- "NA"

#remove all the NA's?
tax16S<- na.omit(tax16S)

#transpose the multinomial, make a column that can be merged by with the taxonomy table
tmultinomial16S <- t(multinomial16S) 
tmultinomial16S <- as.data.frame(tmultinomial16S)
tmultinomial16S <- tibble::rownames_to_column(tmultinomial16S, "row names")
colnames(tmultinomial16S)[which(names(tmultinomial16S) == "row names")] <- "OTUID"

#Do the above with the tax16S table so that there is a column they have in common to merge by
tax16S <- tibble::rownames_to_column(tax16S, "row names")
colnames(tax16S)[which(names(tax16S) == "row names")] <- "OTUID"

#merge the two dataframes 
merged16S_OTUID_Penst <- merge(tmultinomial16S, tax16S, by.x = "OTUID", by.y = "OTUID")
write.csv(merged16S_OTUID_Penst, "merged16S_OTUID_Penst.csv")



##############################################
#### Do the same above with ITS data

metadat <- read.csv("./metadata_wrangled.csv", stringsAsFactors = F)
identity<- metadat[3] #use 22 when you are trying to get the location code in there.
multinomialITS <- read.csv("multinomial_estimates_ITS_by_taxon.csv", stringsAsFactors = F, header = T) # alternatively you can use dirichlet data if you think it's better
multinomialITS<-cbind(identity, multinomialITS)
multinomialITS <- data.frame(multinomialITS[,-1], row.names = multinomialITS[,1])

#standardize by ISD #######WTF is happening with your data
multinomialITS <- multinomialITS[, 2:length(multinomialITS)] /multinomialITS[,names(multinomialITS)== "ISD"]

#remove ISD and duds
multinomialITS<-multinomialITS[,1:(length(multinomialITS)-2)]


#multinomialITS<-multinomialITS[,-2]
taxITS <- import_sintax_file( "penstITS_ESV.sintax", confidence = 0.8)
taxITS<- as.data.frame(taxITS)
taxITS[taxITS == "uncl_Domain"] <- "NA"

#Split the names at the _ so that we can separate out the bacterial names from the syntax deliniations of d,k,p,c,etc.

taxITS$Domain<- gsub("^.*?_", "", taxITS$Domain)
taxITS$Phylum<- gsub("^.*?_", "", taxITS$Phylum)
taxITS$Class<- gsub("^.*?_", "", taxITS$Class)
taxITS$Order<- gsub("^.*?_", "", taxITS$Order)
taxITS$Family<- gsub("^.*?_", "", taxITS$Family)
taxITS$Genus<- gsub("^.*?_", "", taxITS$Genus)
taxITS$Species<- gsub("^.*?_", "", taxITS$Species)

#make sure no unwanted OTU's are present, I did clean my taxonomy but never a bad idea to make sure
taxITS<-subset(taxITS, Family != "mitochondria")
taxITS<-subset(taxITS, Class != "Chloroplast")
taxITS[taxITS == "uncl_Domain"] <- "NA"

#remove all the NA's
taxITS<- na.omit(taxITS)

#transpose data, make a column to merge with the tax table
tmultinomialITS <- t(multinomialITS) # data must be transposed
tmultinomialITS <- as.data.frame(tmultinomialITS)
tmultinomialITS <- tibble::rownames_to_column(tmultinomialITS, "row names")
colnames(tmultinomialITS)[which(names(tmultinomialITS) == "row names")] <- "OTUID"

#Do the above with the tax table so that there is a column they have in common to merge by
taxITS <- tibble::rownames_to_column(taxITS, "row names")
colnames(taxITS)[which(names(taxITS) == "row names")] <- "OTUID"

#merge the dataframes
mergedITS_OTUID_Penst <- merge(tmultinomialITS, taxITS, by.x = "OTUID", by.y = "OTUID")
write.csv(mergedITS_OTUID_Penst, "mergedITS_OTUID_Penst.csv")






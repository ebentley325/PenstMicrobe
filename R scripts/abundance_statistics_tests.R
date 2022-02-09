###################################################
#What Determines Abundance?
###################################################


#Load packages
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scico)
require(RDPutils)

#setwd
setwd("/Users/ebentley/Desktop/PenstMicrobe-master")

#Bring in data--merged16S_OTUID_Penst.csv has everything we need. See the wrangling script or the end of this script for details on formation.
#It has Penstemon sample ID, OTU ID, multinomial data, and OTU taxonomic information
merged16S<- read.csv("./merged16S_OTUID_Penst.csv")
mergedITS<- read.csv("./mergedITS_OTUID_Penst.csv")


#Run lm for each plant taxon
rsquared16S<- list()
k=1
#pdf(file = "16S_abundance~hostTaxa_log10.pdf",width = 10, height = 10)
par(mfrow = c(5,5), mar = c(1,3,1,1), oma = c(2,3,1,1))
for (i in unique(merged16S$Phylum)) {
  #first subset the data
  data_subset <- merged16S[merged16S$Phylum== i, ]
  means <- apply(data_subset[,grep("PE.*", names(data_subset))], 2, FUN=mean)
  host_taxon <- gsub("(\\w+)_\\d+", "\\1", names(means))
  #plot the data in boxplot.log10 to better see the data
  boxplot(log10(means) ~ host_taxon, 
          main = i,
          las = 2,
          xpd = NA, #allow overwriting outside of plot
          xaxt = "n",
          ylab = "abund"
          #ylim =c(0, 6))
  )
  reg <- lm(means ~ host_taxon)
  rsquared16S[[k]]<-summary(reg)$r.squared
  k=k+1
  #text(paste("R^2 ",round(summary(reg)$r.squared,2), sep =""), x =33, y =  par("usr")[4]*0.9, cex=.75)
}
dev.off()

#ITS
rsquaredITS<- list()
k=1
#pdf(file = "ITS_abundance~hostTaxa_log10.pdf",width = 10, height = 10)
par(mfrow = c(3,3), mar = c(1,3,1,1), oma = c(2,3,1,1))
for (i in unique(mergedITS$Phylum)) {
  #first subset the data
  data_subset <- mergedITS[mergedITS$Phylum== i, ]
  means <- apply(data_subset[,grep("PE.*", names(data_subset))], 2, FUN=mean)  
  host_taxon <- gsub("(\\w+)_\\d+", "\\1", names(means))
  #plot
  boxplot(log10(means) ~ host_taxon, 
          main = i,
          las = 2,
          xpd = NA, #allow overwriting outside of plot
          xaxt = "n",
          ylab = "abund"
          #ylim =c(0, 0.15))
  )
  reg <- lm(means ~ host_taxon)
  rsquaredITS[[k]]<-summary(reg)$r.squared
  k=k+1
  #text(paste("R^2 ",round(summary(reg)$r.squared,2), sep =""), x =33, y =  par("usr")[4]*0.9, cex=.75)
}
#dev.off()

















































########### Old code, on here for archival purposes ##########
#### Summarizing Variance in abundances--notes from conversation with Josh ####
#summarize the variance in the abundance of different taxa and ask the question for the taxa that vary a lot:
#--how does host and spatial location relate to abundance (linear models?)
##***## for Chloroflexi 50% of the variance was explained by host taxon--this probably will be afected by read depth. Maybe go back and take only the taxon that have X number of reads
##**## Sequencing depth may be driving the show when it comes to this--JH says almost certainly. 
##**## our options: multinomial helps but the prior will define things and will make them all look the same. The other option is to rarefy and get rid of samples that don't have many reads
##**## It might be worth to do both multinomial data as is and then do it again with rarefied OTU table data with everything under 1000 reads thrown away  
##**## might end up not using the bacterial data bc the read counts are such shit--because of well document sequencing problems having to do with plants


#linear model one phylum at a time as the response variable--a for loop prolly. prolly most interested in beta coefficent for host plant taxa
#whats the effect on this taxon on going from one taxa to another 
#To summarize, could make a heat map or print it out as a mf huge table 
#Maybe all the models will be shit and it won't matter (likely)
#The prior (due to multinomial data) will drive the show and that might be no good at all. Def want to try to repeat the analyses using non-multinomial data--either hellinger transformed data or rarefied data, ISD standardized data, etc

#For each phylum in the data set, run a linear model that uses Taxon and taxa_loc as predictor variables with the abundance as the response variable. In other words,
#For each unique phylum in the dataset, can abundance be predicted by Taxon or taxa_loc?
#i=mergeddata.l[1,]

#calculate the mean abundance of a phylum withing a given penstemon species 

#### Data wrangling 16S Multinomial--can use file "merged16S_OTUID_Penst.csv"--it is what this generates####
#calling both multinomial and taxonomy and doing some formatting for convenient use later
metadat <- read.csv("./metadata_wrangled.csv", stringsAsFactors = F)
identity<- metadat[3] #use 22 when you are trying to get the location code in there.
multinomial <- read.csv("multinomial_estimates_16s_by_taxon.csv", stringsAsFactors = F, header = T) # alternatively you can use dirichlet data if you think it's better
multinomial<-cbind(identity, multinomial)
multinomial <- data.frame(multinomial[,-1], row.names = multinomial[,1])

#standardize by the ISD
multinomial <- multinomial[, 2:length(multinomial)] /multinomial[,(length(multinomial)-1)]

#get rid of the mtDNA, Cholroplast, plant, and duds columns in the multinomial dataframe
multinomial<- multinomial[, 1:(length(multinomial)-4)]

#multinomial<-multinomial[,-3]

tax <- import_sintax_file( "penst16S_ESV.gg.sintax", confidence = 0.8)
tax<- as.data.frame(tax)
tax[tax == "uncl_Domain"] <- "NA"

#Split the names at the _ so that we can separate out the bacterial names from the syntax deliniations of d,k,p,c,etc.

tax$Domain<- gsub("^.*?_", "", tax$Domain)
tax$Phylum<- gsub("^.*?_", "", tax$Phylum)
tax$Class<- gsub("^.*?_", "", tax$Class)
tax$Order<- gsub("^.*?_", "", tax$Order)
tax$Family<- gsub("^.*?_", "", tax$Family)
tax$Genus<- gsub("^.*?_", "", tax$Genus)
tax$Species<- gsub("^.*?_", "", tax$Species)

#make sure no unwanted OTU's are present, I did clean my taxonomy but never a bad idea to make sure
tax<-subset(tax, Family != "mitochondria")
tax<-subset(tax, Class != "Chloroplast")
tax[tax == "uncl_Domain"] <- "NA"

#remove all the NA's?
tax<- na.omit(tax)

#Do not need lines 79-81 when you add in standardizing by the ISD
tmultinomial <- t(multinomial) # data must be transposed
#names <- tmultinomial[1,]
#colnames(tmultinomial) <- names
#tmultinomial <- tmultinomial[-c(1),]

tmultinomial <- as.data.frame(tmultinomial)
tmultinomial <- tibble::rownames_to_column(tmultinomial, "row names")
colnames(tmultinomial)[which(names(tmultinomial) == "row names")] <- "OTUID"

#Do the above with the tax table so that there is a column they have in common to merge by
tax <- tibble::rownames_to_column(tax, "row names")
colnames(tax)[which(names(tax) == "row names")] <- "OTUID"

merged <- merge(tmultinomial, tax, by.x = "OTUID", by.y = "OTUID")
dim(merged)

#### NEW CODE WITH JOSH TO USE NOT LONGFORM DATA ####

#There is variation in bacterial Phylum mean abundance within and among host taxa
#Some plants just have less or more of everything--read counts? 
#16S

rsquared<- list()
k=1
#pdf(file = "16S_for_EPSCoR_Talk.pdf",width = 10, height = 10)
par(mfrow = c(5,5), mar = c(1,3,1,1), oma = c(2,3,1,1))
for (i in unique(merged16S_OTUID_Penst$Phylum)) {
  #first subset the data
  data_subset <- merged16S_OTUID_Penst[merged16S_OTUID_Penst$Phylum== i, ]
  means <- apply(data_subset[,grep("PE.*", names(data_subset))], 2, FUN=median)
  host_taxon <- gsub("(\\w+)_\\d+", "\\1", names(means))
  boxplot(log10(means) ~ host_taxon, 
          main = i,
          las = 2,
          xpd = NA, #allo overwriting outside of plot
          xaxt = "n",
          ylab = "abund"
          #ylim =c(0, 6))
  )
  reg <- lm(means ~ host_taxon)
  rsquared[[k]]<-summary(reg)$r.squared
  k=k+1
  #text(paste("R^2 ",round(summary(reg)$r.squared,2), sep =""), x =33, y =  par("usr")[4]*0.9, cex=.75)
}
#dev.off()

#ITS?
rsquared<- list()
k=1
#pdf(file = "ITS_for_EPSCoR_Talk.pdf",width = 10, height = 10)
par(mfrow = c(3,3), mar = c(1,3,1,1), oma = c(2,3,1,1))
for (i in unique(mergedITS$Phylum)) {
  #first subset the data
  data_subset <- mergedITS[mergedITS$Phylum== i, ]
  means <- apply(data_subset[,grep("PE.*", names(data_subset))], 2, FUN=mean)  
  host_taxon <- gsub("(\\w+)_\\d+", "\\1", names(means))
  boxplot(log10(means) ~ host_taxon, 
          main = i,
          las = 2,
          xpd = NA, #allo overwriting outside of plot
          xaxt = "n",
          ylab = "abund"
          #ylim =c(0, 0.15))
  )
  reg <- lm(means ~ host_taxon)
  rsquared[[k]]<-summary(reg)$r.squared
  k=k+1
  #text(paste("R^2 ",round(summary(reg)$r.squared,2), sep =""), x =33, y =  par("usr")[4]*0.9, cex=.75)
}
#dev.off()




#### older long form code ####

#create a vector that takes sample names
names <- colnames(merged)
tokill <- c("OTUID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
names <- names[!names %in% tokill]
names

mergeddata.l <- reshape2::melt(merged, 
                               id.vars = c("OTUID", "Phylum", "Class", "Order", "Family", "Genus"), 
                               measure.vars = names ,
                               value.name = "abundance", 
                               variable.name = "site")

mergeddata.l$abundance <- as.numeric(mergeddata.l$abundance)
#Save this table so that you can use it without having to run through all this code every time


#### 16S multinomial plots of abundance and R2 ####
#First, use gsub to make a plant taxon field
#mergeddata.l$host_taxon<-gsub("($.*_)", "\\1", mergeddata.l$site)


mergeddata.l$host_taxon <- gsub("(\\w+)_\\d+", "\\1", mergeddata.l$site)


rsquared<- list()
k=1
par(mfrow = c(5,5), mar = c(1,3,1,1), oma = c(2,3,1,1))
for (i in unique(mergeddata.l$Phylum)) {
  #first subset the data
  data_subset <- mergeddata.l[mergeddata.l$Phylum== i, ]
  boxplot(log10(data_subset$abundance) ~ data_subset$host_taxon, 
          main = i,
          las = 2,
          xpd = NA, #allo overwriting outside of plot
          xaxt = "n",
          ylab = "abund")
  reg <- lm(data_subset$abundance ~ data_subset$host_taxon)
  rsquared[[k]]<-summary(reg)$r.squared
  k=k+1
  #text(paste("R^2 ",round(summary(reg)$r.squared,2), sep =""), x =33, y =  par("usr")[4]*0.9, cex=.75)
}



penstemon_phylum_medianAbund <- aggregate(mergeddata.l$abundance ~ mergeddata.l$host_taxon + mergeddata.l$Phylum, FUN = median)

options(scipen = 99)

#Determine the most abundant phyla
penstemon_phylum_mostabund <- aggregate(penstemon_phylum_medianAbund$`mergeddata.l$abundance` ~ penstemon_phylum_medianAbund$`mergeddata.l$Phylum`, FUN = median)
penstemon_phylum_mostabund[ order(penstemon_phylum_mostabund[,2]), ]


#### ITS Multinomial data wrangling ####

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

#remove all the NA's?
taxITS<- na.omit(taxITS)

#Don't need names/colnames/tmulti if standardized by ISD
tmultinomialITS <- t(multinomialITS) # data must be transposed
#names <- tmultinomialITS[1,]
#colnames(tmultinomialITS) <- names
#tmultinomialITS <- tmultinomialITS[-c(1),]

tmultinomialITS <- as.data.frame(tmultinomialITS)
tmultinomialITS <- tibble::rownames_to_column(tmultinomialITS, "row names")
colnames(tmultinomialITS)[which(names(tmultinomialITS) == "row names")] <- "OTUID"

#Do the above with the tax table so that there is a column they have in common to merge by
taxITS <- tibble::rownames_to_column(taxITS, "row names")
colnames(taxITS)[which(names(taxITS) == "row names")] <- "OTUID"

mergedITS <- merge(tmultinomialITS, taxITS, by.x = "OTUID", by.y = "OTUID")

#create a vector that takes sample names
names <- colnames(mergedITS)
tokill <- c("OTUID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
names <- names[!names %in% tokill]
names

mergeddata.l.ITS <- reshape2::melt(mergedITS, 
                               id.vars = c("OTUID", "Phylum", "Class", "Order", "Family", "Genus"), 
                               measure.vars = names ,
                               value.name = "abundance", 
                               variable.name = "site")

mergeddata.l.ITS$abundance <- as.numeric(mergeddata.l.ITS$abundance)


#### ITS multinomial plots of abundance and R2 ####
#First, use gsub to make a plant taxon field
#mergeddata.l.ITS$host_taxon<-gsub("($.*_)", "\\1", mergeddata.l.ITS$site)


mergeddata.l.ITS$host_taxon <- gsub("(\\w+)_\\d+", "\\1", mergeddata.l.ITS$site)


rsquaredITS<- list()
k=1
pdf(file = "ITS_abundance~HostTaxon_ISD_Standardized.pdf",width = 10, height = 10)
par(mfrow = c(3,3), mar = c(1,3,1,1), oma = c(2,4,1,1))
for (i in unique(mergeddata.l.ITS$Phylum)) {
  #first subset the data
  data_subset <- mergeddata.l.ITS[mergeddata.l.ITS$Phylum== i, ]
 boxplot(log10(data_subset$abundance) ~ data_subset$host_taxon, 
          main = i,
          las = 2,
          xpd = NA, #allow overwriting outside of plot
          xaxt = "n",
          ylab = "abund")
 
 reg <- lm(data_subset$abundance ~ data_subset$host_taxon)
  rsquaredITS[[k]]<-summary(reg)$r.squared
  k=k+1
  #plot the R2 on the plots
  #text(paste("R^2 ",round(summary(reg)$r.squared,2), sep =""), x =30, y =  par("usr")[4]*0.9, cex=1)
}
dev.off()



penstemon_phylum_medianAbund <- aggregate(mergeddata.l$abundance ~ mergeddata.l$host_taxon + mergeddata.l$Phylum, FUN = median)

options(scipen = 99)

#Determine the most abundant phyla
penstemon_phylum_mostabund <- aggregate(penstemon_phylum_medianAbund$`mergeddata.l$abundance` ~ penstemon_phylum_medianAbund$`mergeddata.l$Phylum`, FUN = median)
penstemon_phylum_mostabund[ order(penstemon_phylum_mostabund[,2]), ]


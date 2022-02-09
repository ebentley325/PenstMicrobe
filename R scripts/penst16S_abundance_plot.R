#### USING SEIF'S CODE THAT WILL HOPEFULLY WORK BETTER ####

#Load packages
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scico)
require(RDPutils)

#setwd
setwd("/Users/ebentley/Desktop/PenstMicrobe-master")

#calling both multinomial and taxonomy and doing some formatting for convenient use later
metadat <- read.csv("./metadata_wrangled.csv", stringsAsFactors = F)
identity<- metadat[,c(3,22)]
multinomial <- read.csv("multinomial_estimates_16s_by_taxon.csv", stringsAsFactors = F, header = T) # alternatively you can use dirichlet data if you think it's better
multinomial<-cbind(identity, multinomial)
#get rid of the mtDNA, Cholroplast, plant, and duds columns in the multinomial dataframe
multinomial<- multinomial[, 1:(length(multinomial)-4)]

#multinomial <- data.frame(multinomial[,-1], row.names = multinomial[,1])
multinomial<-multinomial[,-3]
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

tmultinomial <- t(multinomial) # data must be transposed
names <- tmultinomial[1,]
colnames(tmultinomial) <- names
tmultinomial <- tmultinomial[-c(1),]

tmultinomial <- as.data.frame(tmultinomial)
tmultinomial <- tibble::rownames_to_column(tmultinomial, "row names")
colnames(tmultinomial)[which(names(tmultinomial) == "row names")] <- "OTUID"

#Do the above with the tax table so that there is a column they have in common to merge by
tax <- tibble::rownames_to_column(tax, "row names")
colnames(tax)[which(names(tax) == "row names")] <- "OTUID"

merged <- merge(tmultinomial, tax, by.x = "OTUID", by.y = "OTUID")
dim(merged)


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

#### SEIF PLOT! PLOT! PLOT! ####

#the most basic non polished crappy looking graph
ggplot(mergeddata.l, aes(fill=Family, y=abundance, x=site)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_discrete(guide = "none")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

#Another, different graph
ggplot(mergeddata.l, aes(x = site, y = abundance)) +
  facet_wrap(~Phylum, scales = "free") +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

#This next boy did NOT want to work with me
#selecting few 
ggplot(subset(mergeddata.l, Phylum %in% c("Acidobacteria",  "Actinobacteria")), #, "Bacteroidetes", "Chloroflexi"
       aes(x = site, y = abundance)) +
  facet_wrap(~Genus, scales = "free") +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

#How many unique Phyla do we have
phylaUnique<-unique(mergeddata.l[c("Phylum")])
#24 unique phyla

#How many unique genera do we have
generaUnique<-unique(mergeddata.l[c("Genus")])
#278 unique genera

#### Summarize data so we can polish figures ####

## create column with sample name (= site name without the technical replicate id)
#this doesn't actually do anything for me but could do something similar to separate taxa and loc
#mergeddata.l$sample <- mergeddata.l$site
#mergeddata.l$sample <- gsub("\\.\\d+", "", mergeddata.l$sample)


## calculate mean per sample (ie mean of technical replicates): this is to get rid of technical replicates, again you might have already done this for your data, I did means of technical replicates, some other folks did sum of both
#we merged technical replicates BEFORE CNVRG modelling
#mergeddata.l.mean <- mergeddata.l %>% 
 # group_by(OTUID, Phylum, Class, Order, Family, Genus, sample) %>%
 # summarise(abundance_mean = mean(abundance, na.rm = TRUE))

#wide format of dataframe (OTU like table): this was for me to save table with merged technical replicates
mergeddata.l.mean.wide <- dcast(mergeddata.l.mean, sample~OTUID, value.var = "abundance_mean")

pdf(width = 7, height = 7, file = "abundance_mean_all_samples.pdf")
ggplot(mergeddata.l.mean, aes(x = sample, y = abundance_mean)) +
  facet_wrap(~Phylum, scales = "free") +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
dev.off()

pdf(width = 15, height = 7, file = "abundance_stacked_barplot_all_samples")
ggplot(mergeddata.l.mean, aes(fill=Phylum, y=abundance_mean, x=sample)) + 
  geom_bar(position="fill", stat="identity")+
  # scale_fill_discrete(guide = "none")+ # this or viridis
  scale_fill_viridis_d(option = "inferno")+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=5)
  )
dev.off()

# 10 most abundant phyla, now I want my figure to reflect only top 10 phyla --------------------------------------------------
# sum phylum abundance by sample
# If you want a figure with other taxonomic class than phylum, simply replace in group_by below
summed_phyla_sample <-mergeddata.l.mean %>% 
  group_by(sample, Phylum) %>%
  summarise(abundance_sum = sum(abundance_mean, na.rm = T)) %>% 
  arrange(desc(abundance_sum))

summed_phyla_sample<- na.omit(summed_phyla_sample)

sample_phyla_abundance_table<-as.data.frame(summed_phyla_sample)
write.csv(sample_phyla_abundance_table, "sample_phyla_abundance_table")

top10 <- summed_phyla_sample[1:10,]

top10

ggplot(subset(top10, sample!="Blank"), aes(fill=Phylum, y=abundance_sum, x=sample)) + 
  geom_bar(position="fill", stat="identity", colour = "lightgray")+
  #scale_fill_scico_d(palette = "oleron") +
  scale_fill_brewer(palette = "Paired")+ # this or viridis
  #scale_fill_viridis_d(option = "viridis")+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)
  ) +
  ggtitle("ten most abundant prokaryotic phyla per sample") +
  ylab("Proportional abundance") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())

ggsave("prokaryotes.png", width = 32, height = 15, units = "cm")



#### P. Haydenii ####
PHay_multi<-multinomial[grep("PEHA12", multinomial$identity, perl=TRUE), ]

PHay_tmultinomial <- t(PHay_multi) # data must be transposed
names <- PHay_tmultinomial[1,]
colnames(PHay_tmultinomial) <- names
PHay_tmultinomial <- PHay_tmultinomial[-c(1),]

PHay_tmultinomial <- as.data.frame(PHay_tmultinomial)
PHay_tmultinomial <- tibble::rownames_to_column(PHay_tmultinomial, "row names")
colnames(PHay_tmultinomial)[which(names(PHay_tmultinomial) == "row names")] <- "OTUID"

#Do the above with the tax table so that there is a column they have in common to merge by
#only do this once
#tax <- tibble::rownames_to_column(tax, "row names")
#colnames(tax)[which(names(tax) == "row names")] <- "OTUID"

PHay_merged <- merge(PHay_tmultinomial, tax, by.x = "OTUID", by.y = "OTUID")
dim(PHay_merged)


#create a vector that takes sample names
names <- colnames(PHay_merged)
tokill <- c("OTUID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
names <- names[!names %in% tokill]
names

Phay_mergeddata.l <- reshape2::melt(PHay_merged, 
                               id.vars = c("OTUID", "Phylum", "Class", "Order", "Family", "Genus"), 
                               measure.vars = names ,
                               value.name = "abundance", 
                               variable.name = "site")

Phay_mergeddata.l$abundance <- as.numeric(Phay_mergeddata.l$abundance)

Phay_mergeddata.l$sample <- Phay_mergeddata.l$site
Phay_mergeddata.l$sample <- gsub("\\.\\d+", "", Phay_mergeddata.l$sample)


## calculate mean per sample (ie mean of technical replicates): this is to get rid of technical replicates, again you might have already done this for your data, I did means of technical replicates, some other folks did sum of both
PHay_mergeddata.l.mean <- Phay_mergeddata.l %>% 
  group_by(OTUID, Phylum, Class, Order, Family, Genus, sample) %>%
  summarise(abundance_mean = mean(abundance, na.rm = TRUE))

#wide format of dataframe (OTU like table): this was for me to save table with merged technical replicates
#mergeddata.l.mean.wide <- dcast(PHay_mergeddata.l.mean, sample~OTUID, value.var = "abundance_mean")

pdf(width = 7, height = 7, file = "PEHA12_stacked_abundance_plot.pdf")
ggplot(PHay_mergeddata.l.mean, aes(fill=Phylum, y=abundance_mean, x=sample)) + 
  geom_bar(position="fill", stat="identity")+
  # scale_fill_discrete(guide = "none")+ # this or viridis
  scale_fill_viridis_d(option = "inferno")+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=5)
  )
dev.off()


#### PEGR7 ####
PGrand_multi<-multinomial[grep("PEGR7", multinomial$identity, perl=TRUE), ]

PGrand_tmultinomial <- t(PGrand_multi) # data must be transposed
names <- PGrand_tmultinomial[1,]
colnames(PGrand_tmultinomial) <- names
PGrand_tmultinomial <- PGrand_tmultinomial[-c(1),]

PGrand_tmultinomial <- as.data.frame(PGrand_tmultinomial)
PGrand_tmultinomial <- tibble::rownames_to_column(PGrand_tmultinomial, "row names")
colnames(PGrand_tmultinomial)[which(names(PGrand_tmultinomial) == "row names")] <- "OTUID"

#Do the above with the tax table so that there is a column they have in common to merge by. Only do it here if you 
#haven't run the code above ERIN
#tax <- tibble::rownames_to_column(tax, "row names")
#colnames(tax)[which(names(tax) == "row names")] <- "OTUID"

PGrand_merged <- merge(PGrand_tmultinomial, tax, by.x = "OTUID", by.y = "OTUID")
#dim(merged)


#create a vector that takes sample names
names <- colnames(PGrand_merged)
tokill <- c("OTUID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
names <- names[!names %in% tokill]
names

PGrand_mergeddata.l <- reshape2::melt(PGrand_merged, 
                               id.vars = c("OTUID", "Phylum", "Class", "Order", "Family", "Genus"), 
                               measure.vars = names ,
                               value.name = "abundance", 
                               variable.name = "site")

PGrand_mergeddata.l$abundance <- as.numeric(PGrand_mergeddata.l$abundance)

PGrand_mergeddata.l$sample <- PGrand_mergeddata.l$site
PGrand_mergeddata.l$sample <- gsub("\\.\\d+", "", PGrand_mergeddata.l$sample)


## calculate mean per sample (ie mean of technical replicates): this is to get rid of technical replicates, again you might have already done this for your data, I did means of technical replicates, some other folks did sum of both
PGrand_mergeddata.l.mean <- PGrand_mergeddata.l %>% 
  group_by(OTUID, Phylum, Class, Order, Family, Genus, sample) %>%
  summarise(abundance_mean = mean(abundance, na.rm = TRUE))

#wide format of dataframe (OTU like table): this was for me to save table with merged technical replicates
#PGrand_mergeddata.l.mean <- dcast(PGrand_mergeddata.l.mean, sample~OTUID, value.var = "abundance_mean")

pdf(width = 7, height = 7, file = "PEGR7_stacked_abundance_plot.pdf")
ggplot(PGrand_mergeddata.l.mean, aes(fill=Phylum, y=abundance_mean, x=sample)) + 
  geom_bar(position="fill", stat="identity")+
  # scale_fill_discrete(guide = "none")+ # this or viridis
  scale_fill_viridis_d(option = "inferno")+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
  )
dev.off()

#### PECA 17 ####
PECA17_multi<-multinomial[grep("PECA", multinomial$identity, perl=TRUE), ]

PECA17_tmulti <- t(PECA17_multi) # data must be transposed
names <- PECA17_tmulti[1,]
colnames(PECA17_tmulti) <- names
PECA17_tmulti <- PECA17_tmulti[-c(1),]

PECA17_tmulti <- as.data.frame(PECA17_tmulti)
PECA17_tmulti <- tibble::rownames_to_column(PECA17_tmulti, "row names")
colnames(PECA17_tmulti)[which(names(PECA17_tmulti) == "row names")] <- "OTUID"

#Do the above with the tax table so that there is a column they have in common to merge by. Only do it here if you 
#haven't run the code above ERIN
#tax <- tibble::rownames_to_column(tax, "row names")
#colnames(tax)[which(names(tax) == "row names")] <- "OTUID"

PECA17_merged <- merge(PECA17_tmulti, tax, by.x = "OTUID", by.y = "OTUID")
#dim(merged)


#create a vector that takes sample names
names <- colnames(PECA17_merged)
tokill <- c("OTUID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
names <- names[!names %in% tokill]
names

PECA17_mergeddata.l <- reshape2::melt(PECA17_merged, 
                                      id.vars = c("OTUID", "Phylum", "Class", "Order", "Family", "Genus"), 
                                      measure.vars = names ,
                                      value.name = "abundance", 
                                      variable.name = "site")

PECA17_mergeddata.l$abundance <- as.numeric(PECA17_mergeddata.l$abundance)

PECA17_mergeddata.l$sample <- PECA17_mergeddata.l$site
PECA17_mergeddata.l$sample <- gsub("\\.\\d+", "", PECA17_mergeddata.l$sample)


## calculate mean per sample (ie mean of technical replicates): this is to get rid of technical replicates, again you might have already done this for your data, I did means of technical replicates, some other folks did sum of both
PECA_mergeddata.l.mean <- PECA17_mergeddata.l %>% 
  group_by(OTUID, Phylum, Class, Order, Family, Genus, sample) %>%
  summarise(abundance_mean = mean(abundance, na.rm = TRUE))

#wide format of dataframe (OTU like table): this was for me to save table with merged technical replicates
#PGrand_mergeddata.l.mean <- dcast(PGrand_mergeddata.l.mean, sample~OTUID, value.var = "abundance_mean")

pdf(width = 7, height = 7, file = "PECA17_stacked_abundance_plot.pdf")
ggplot(PECA_mergeddata.l.mean, aes(fill=Phylum, y=abundance_mean, x=sample)) + 
  geom_bar(position="fill", stat="identity")+
  # scale_fill_discrete(guide = "none")+ # this or viridis
  scale_fill_viridis_d(option = "inferno")+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
        )
dev.off()

# PECA17_49 is weeeeird 
# Does it have a weird number of reads?

OTU_table<- read.csv("./16S_otuTable_arranged_by_taxon.csv")
OTU_table<-OTU_table[,1:1430]
OTU_table$sums<- rowSums(OTU_table[,2:1430])

#Make the first column the rownames
OTU_table <- data.frame(OTU_table[,-1], row.names = OTU_table[,1])
OTU_table <- tibble::rownames_to_column(OTU_table, "row names")
colnames(OTU_table)[which(names(OTU_table) == "row names")] <- "sample"

weird <- OTU_table[OTU_table$sample == "PECA17_49",]
sum(weird[, 2:(length(weird)-5)])
table(weird[, 2:(length(weird)-5)] > 1)

pdf(width = 7, height = 7, file = "readcount16S.pdf")
ggplot(OTU_table, aes(y=sums, x=sample)) + 
  geom_col() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7) 
  )
dev.off()


## HOW DOES ITS LOOK
ITS_Taxon_OTU<- read.csv("./ITS_otuTable_arranged_by_taxon.csv")
ITS_Taxon_OTU<-ITS_Taxon_OTU[,1:7238]
ITS_Taxon_OTU$sums<- rowSums(ITS_Taxon_OTU[,2:7238])

#Make the first column the rownames
ITS_Taxon_OTU <- data.frame(ITS_Taxon_OTU[,-1], row.names = ITS_Taxon_OTU[,1])
ITS_Taxon_OTU <- tibble::rownames_to_column(ITS_Taxon_OTU, "row names")
colnames(ITS_Taxon_OTU)[which(names(ITS_Taxon_OTU) == "row names")] <- "sample"

pdf(width = 7, height = 7, file = "readCountITS.pdf")
ggplot(ITS_Taxon_OTU, aes(y=sums, x=sample)) + 
  geom_col() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7) 
  )
dev.off()

#########################################################################################################################################################################

#### Abundance plots for ITS ####

multinomial <- read.csv("multinomial_estimates_ITS_by_taxon.csv", stringsAsFactors = F, header = T) # alternatively you can use dirichlet data if you think it's better
multinomial<-cbind(identity, multinomial)
multinomial <- data.frame(multinomial[,-1], row.names = multinomial[,1])
multinomial<-multinomial[,-2]
tax <- import_sintax_file( "penstITS_ESV.sintax", confidence = 0.8)
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

tmultinomial <- t(multinomial) # data must be transposed
names <- tmultinomial[1,]
colnames(tmultinomial) <- names
tmultinomial <- tmultinomial[-c(1),]

tmultinomial <- as.data.frame(tmultinomial)
tmultinomial <- tibble::rownames_to_column(tmultinomial, "row names")
colnames(tmultinomial)[which(names(tmultinomial) == "row names")] <- "OTUID"

#Do the above with the tax table so that there is a column they have in common to merge by
tax <- tibble::rownames_to_column(tax, "row names")
colnames(tax)[which(names(tax) == "row names")] <- "OTUID"

merged <- merge(tmultinomial, tax, by.x = "OTUID", by.y = "OTUID")
dim(merged)


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

## create column with sample name (= site name without the technical replicate id)
mergeddata.l$sample <- mergeddata.l$site
mergeddata.l$sample <- gsub("\\.\\d+", "", mergeddata.l$sample)


## calculate mean per sample (ie mean of technical replicates): this is to get rid of technical replicates, again you might have already done this for your data, I did means of technical replicates, some other folks did sum of both
mergeddata.l.mean <- mergeddata.l %>% 
  group_by(OTUID, Phylum, Class, Order, Family, Genus, sample) %>%
  summarise(abundance_mean = mean(abundance, na.rm = TRUE))

pdf(width = 12, height = 7, file = "abundance_stacked_barplot_ITS.pdf")
ggplot(mergeddata.l.mean, aes(fill=Phylum, y=abundance_mean, x=sample)) + 
  geom_bar(position="fill", stat="identity")+
  # scale_fill_discrete(guide = "none")+ # this or viridis
  scale_fill_viridis_d(option = "inferno")+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=5)
  )
dev.off()




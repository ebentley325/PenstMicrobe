rm(list=ls())
multi <- read.csv("multinomial_estimates_ITS_by_taxon.csv", 
                  stringsAsFactors = F)
dim(multi)
metadat <- read.csv("./metadata_wrangled.csv", stringsAsFactors = F)
multi[,1] <- metadat$Identifier

#Looks like blanks were removed already. 

capscale()
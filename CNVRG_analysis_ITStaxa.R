rm(list=ls())
library(CNVRG)
dat <- read.csv("./ITS_otuTable_arranged_by_taxon.csv",
                  fill = T, header = T, stringsAsFactors = F)

metadat <- read.csv("./metadata_wrangled.csv", stringsAsFactors = F)
table(dat$X == metadat$Identifier)
table(metadat$Taxon)

#add a one to our data
#also remove any zeros
#no zeros
table(colSums(dat[,2:(length(dat)-2)]) == 0)
dat[,2:length(dat)] <- dat[,2:length(dat)] + 1

#Metal AF!
set.seed(666)

modelOut <- cnvrg_VI(
  countData = dat,
  starts = indexer(metadat$Taxon)$starts,
  ends = indexer(metadat$Taxon)$ends,
  output_samples = 500,
  params_to_save = c("pi", "p")
)
save.image(file = "./CNVRG_ITStaxon.Rdata")
diffs <- diff_abund(model_output = modelOut, countData = dat)

save.image(file = "./CNVRG_ITSdiffstaxon.Rdata")


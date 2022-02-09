#### PCoA for ITS Values ####

rm(list=ls())
library(RColorBrewer)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

multi <- read.csv("multinomial_estimates_ITS_by_taxon.csv", 
                  stringsAsFactors = F)
dim(multi)
metadat <- read.csv("./metadata_wrangled.csv", stringsAsFactors = F)
multi[,1] <- metadat$Identifier

metadat$colors <- NA

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, 
                            qual_col_pals$maxcolors, 
                            rownames(qual_col_pals)))

k <- 1
for(i in unique(metadat$Taxon)){
  metadat$colors[metadat$Taxon == i] <- col_vector[k]
  k <- k + 1
}
metadat$colors[metadat$colors == "#FFFF99"] <- "darkorchid4"
metadat$colors[metadat$colors == "#B2DF8A"] <- "brown"

metadat$shapes <- 21
#metadat$shapes[metadat$colors == "#B2DF8A"] <- 24
#metadat$shapes[metadat$colors == "#B3E2CD"] <- 8
#metadat$shapes[metadat$colors == "#CCEBC5"] <- 11

#blanks were removed already. 

#CHANGE THE index in the divisor so that it is the ISD
#this makes it standardized by the ISD. Need to subtract more for the 16S datasets

multi <- multi[, 2:length(multi)] /multi[,(length(multi)-1)]
multi<- multi[,1:376]
#Change the indexing so that we remove, plant, mtdna, isd, duds. Probably
#wll need to subtract 4, but not sure
#EB--only -2, the last two columns are "ISD" and "Duds"
# In the multi16S data the last four will need to be removed as there are columns
#labeled plant, IDS, duds, mtDNA

bcdist <- ecodist::bcdist(multi[, 2:(length(multi)-2)]) #Bray-Curtis dissimilarity

PCOA <- cmdscale(d = bcdist, k = 4) #Do the analysis

#pdf(width = 7, height = 7, file = "pcoa_by_taxon_div_isd.pdf")
plot(type = "n",PCOA[,1], PCOA[,2],
     ylab = "PCoA #2 (15%)",
     frame.plot = F,
    # ylim = c(-.4,0.4),
    # xlim = c(-.4,0.4),
     xlab = "PCoA #1 (23%)",
     yaxt = "n",
     xaxt = "n",
     main = "ITS by Taxon",
     las = 2)
axis(side = 1, at = seq(-0.4,0.4,0.2), lab = seq(-0.4,0.4,0.2))
axis(side = 2, at = seq(-0.4,0.4,0.2), lab = seq(-0.4,0.4,0.2), las = 2)

points(PCOA[,1], PCOA[,2],
       pch = metadat$shapes,
       cex = 1.8,
       xpd = NA,
       bg = add.alpha(metadat$colors, 0.7),
       col = add.alpha(metadat$colors, 0.9))
#dev.off()

pdf(width = 8, height = 5, file = "pcoa_by_taxon_legend_current.pdf")
#Legend plot
plot(type = "n",
     PCOA[,1], PCOA[,2],
    # ylim = c(0,10),
     axes = F,
     frame.plot = F,
     xlab = "",
     ylab = ""
    #xlim = c(0,4)
     )
legend(x = 0, y = 10, 
       #legend = unique(metadat$colors),
       legend = unique(metadat$Taxon),
       cex = 0.8,
       pt.cex = 1.4,
       xpd = NA,
       bty = "n",
       y.intersp=1,
       ncol = 2,
       pch =  16, 
       col = unique(metadat$colors))
dev.off()

#extracting the variance explained by each axis
PCOAeigs <- cmdscale(d = bcdist, k = 4, eig = TRUE) #Do the analysis

#variance for first axis
eigs <- PCOAeigs$eig
eigs[1] / sum(eigs)

#variance for second
eigs[2] / sum(eigs)







#### PCoA for 16S ####
#make 16S into a readable format
multi16S <- read.csv("multinomial_estimates_16s_by_taxon.csv", 
                  stringsAsFactors = F)
dim(multi16S)

#replace with taxon identifiers
multi16S[,1] <- metadat$Identifier

#Make it so that the index is the taxon identifier
rownames(multi16S) <- multi16S$X
multi16S<- multi16S[,2:1434]


#CHANGE THE index in the divisor so that it is the ISD

multi16S <- multi16S[, 2:length(multi16S)] /multi16S[,(length(multi16S)-1)]

#Remove 4 to get rid of plant, ISD, duds, and mtDNA

bcdist16S <- ecodist::bcdist(multi16S[, 2:(length(multi16S)-4)]) #Bray-Curtis dissimilarity

PCOA16S <- cmdscale(d = bcdist16S, k = 4) #Do the analysis

#pdf(width = 7, height = 7, file = "pcoa16S_by_taxon_div_isd.pdf")
plot(type = "n",PCOA16S[,1],PCOA16S[,2],
    # ylab = "PCoA #2 (20%)",
     frame.plot = F,
     ylim = c(-.6,0.4),
     xlim = c(-.4,0.6),
     #xlab = "PCoA #1 (58%)",
     yaxt = "n",
     xaxt = "n",
     main = "16S by Taxon",
     las = 2)
axis(side = 1, at = seq(-0.4,0.6,0.2), lab = seq(-0.4,0.6,0.2))
axis(side = 2, at = seq(-0.6,0.4,0.2), lab = seq(-0.6,0.4,0.2), las = 2)

###*****fix this graph so that that one point isn't below the axis
points(PCOA16S[,1], PCOA16S[,2],
       pch = metadat$shapes,
       cex = 1.8,
       xpd = NA,
       bg = add.alpha(metadat$colors, 0.7),
       col = add.alpha(metadat$colors, 0.9))
#dev.off()


pdf(width = 8, height = 5, file = "pcoa16S_by_taxon_legend_current.pdf")
#Legend plot
plot(type = "n",
     PCOA16S[,1], PCOA16S[,2],
     # ylim = c(0,10),
     axes = F,
     frame.plot = F,
     xlab = "",
     ylab = ""
     #xlim = c(0,4)
)
legend(x = 0, y = 10, 
       #legend = unique(metadat$colors),
       legend = unique(metadat$Taxon),
       cex = 0.8,
       pt.cex = 1.4,
       xpd = NA,
       bty = "n",
       y.intersp=1,
       ncol = 2,
       pch =  16, 
       col = unique(metadat$colors))
dev.off()

#Variance explained calculation

PCOA16Seigs <- cmdscale(d = bcdist16S, k = 4, eig = TRUE) #Do the analysis

#variance of axis 1
eigs16S <- PCOA16Seigs$eig
eigs16S[1] / sum(eigs16S)

#variance of axis 2
eigs16S[2] / sum(eigs16S)


#### PCOA for taxon + location ####
multiITStax_site <- read.csv("multinomial_estimates_ITS_by_taxon_site.csv", 
                  stringsAsFactors = F)

multiITStax_site[,1] <- metadat$Identifier

multiITStax_site <- multiITStax_site[, 2:length(multiITStax_site)] /multiITStax_site[,(length(multiITStax_site)-1)]

#Remove 2 to get rid of IDS and duds
bcdist_tax_site <- ecodist::bcdist(multiITStax_site[, 2:(length(multiITStax_site)-2)]) #Bray-Curtis dissimilarity

PCOAITStax_site <- cmdscale(d = bcdist_tax_site, k = 4) #Do the analysis
#pdf(width = 7, height = 7, file = "pcoaITS_by_taxon_site_div_isd.pdf")
plot(type = "n",PCOAITStax_site[,1],PCOAITStax_site[,2],
     ylab = "PCoA #2 (15%)",
     frame.plot = F,
     ylim = c(-.3,0.5),
     xlim = c(-.5,0.5),
     xlab = "PCoA #1 (22%)",
     yaxt = "n",
     xaxt = "n",
     main = "ITS by Taxon and Site",
     las = 2)
axis(side = 1, at = seq(-0.5,0.5,0.2), lab = seq(-0.5,0.5,0.2))
axis(side = 2, at = seq(-0.3,0.5,0.2), lab = seq(-0.3,0.5,0.2), las = 2)

points(PCOAITStax_site[,1], PCOAITStax_site[,4],
       pch = metadat$shapes,
       cex = 1.8,
       xpd = NA,
       bg = add.alpha(metadat$colors, 0.7),
       col = add.alpha(metadat$colors, 0.9))
#dev.off()


pdf(width = 8, height = 5, file = "pcoaITS_by_taxon_site_legend_current.pdf")
#Legend plot
plot(type = "n",
     PCOAITStax_site[,1], PCOAITStax_site[,2],
     # ylim = c(0,10),
     axes = F,
     frame.plot = F,
     xlab = "",
     ylab = ""
     #xlim = c(0,4)
)
legend(x = 0, y = 10, 
       #legend = unique(metadat$colors),
       legend = unique(metadat$Taxon),
       cex = 0.8,
       pt.cex = 1.4,
       xpd = NA,
       bty = "n",
       y.intersp=1,
       ncol = 2,
       pch =  16, 
       col = unique(metadat$colors))
dev.off()

#calculate variance
PCOAITStax_site_eigs <- cmdscale(d = bcdist_tax_site, k = 4, eig = TRUE) #Do the analysis

#Axis 1
eigITS_tax_site <- PCOAITStax_site_eigs$eig
eigITS_tax_site[1] / sum(eigITS_tax_site)

#axis 2
eigITS_tax_site[2] / sum(eigITS_tax_site)



#### 16s by Taxon and site ####
multi16Stax_site <- read.csv("multinomial_estimates_16s_by_taxon_site.csv", 
                             stringsAsFactors = F)

multi16Stax_site[,1] <- metadat$Identifier

multi16Stax_site <- multi16Stax_site[, 2:length(multi16Stax_site)] /multi16Stax_site[,(length(multi16Stax_site)-1)]

#Remove 4
bcdist_16Stax_site <- ecodist::bcdist(multi16Stax_site[, 2:(length(multi16Stax_site)-4)]) #Bray-Curtis dissimilarity

PCOA16Stax_site <- cmdscale(d = bcdist_16Stax_site, k = 4) #Do the analysis
#pdf(width = 7, height = 7, file = "pcoaITS_by_taxon_site_div_isd.pdf")
plot(type = "n",PCOA16Stax_site[,1],PCOA16Stax_site[,2],
     ylab = "PCoA #2 (15%)",
     frame.plot = F,
     ylim = c(-.2,0.6),
     xlim = c(-.4,0.6),
     xlab = "PCoA #1 (22%)",
     yaxt = "n",
     xaxt = "n",
     main = "16S by Taxon and Site",
     las = 2)
axis(side = 1, at = seq(-0.4,0.6,0.2), lab = seq(-0.4,0.6,0.2))
axis(side = 2, at = seq(-0.2,0.6,0.2), lab = seq(-0.2,0.6,0.2), las = 2)

points(PCOA16Stax_site[,1], PCOA16Stax_site[,4],
       pch = metadat$shapes,
       cex = 1.8,
       xpd = NA,
       bg = add.alpha(metadat$colors, 0.7),
       col = add.alpha(metadat$colors, 0.9))
#dev.off()


pdf(width = 8, height = 5, file = "pcoa16S_by_taxon_site_legend_current.pdf")
#Legend plot
plot(type = "n",
     PCOA16Stax_site[,1], PCOA16Stax_site[,2],
     # ylim = c(0,10),
     axes = F,
     frame.plot = F,
     xlab = "",
     ylab = ""
     #xlim = c(0,4)
)
legend(x = 0, y = 10, 
       #legend = unique(metadat$colors),
       legend = unique(metadat$taxa_loc),
       cex = 0.8,
       pt.cex = 1.4,
       xpd = NA,
       bty = "n",
       y.intersp=1,
       ncol = 2,
       pch =  16, 
       col = unique(metadat$colors))
dev.off()
#variance explained

PCOAITStax_site_eigs <- cmdscale(d = bcdist_tax_site, k = 4, eig = TRUE) #Do the analysis

#Axis 1
eigITS_tax_site <- PCOAITStax_site_eigs$eig
eigITS_tax_site[1] / sum(eigITS_tax_site)

#axis 2
eigITS_tax_site[2] / sum(eigITS_tax_site)



#### MAKE THE LEGEND PLOTS WORK





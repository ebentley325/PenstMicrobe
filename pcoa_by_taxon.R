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
metadat$shapes[metadat$colors == "#B2DF8A"] <- 24
metadat$shapes[metadat$colors == "#B3E2CD"] <- 8
metadat$shapes[metadat$colors == "#CCEBC5"] <- 11

#Looks like blanks were removed already. 
bcdist <- ecodist::bcdist(multi[,2:(length(multi)-2)]) #Bray-Curtis dissimilarity

PCOA <- cmdscale(d = bcdist, k = 4) #Do the analysis

pdf(width = 7, height = 7, file = "pcoa_by_taxon.pdf")
plot(type = "n",PCOA[,1], PCOA[,2],
     ylab = "PCoA #2",
     frame.plot = F,
     ylim = c(-.4,0.4),
     xlim = c(-.4,0.4),
     xlab = "PCoA #1",
     yaxt = "n",
     xaxt = "n",
     las = 2)
axis(side = 1, at = seq(-0.4,0.4,0.2), lab = seq(-0.4,0.4,0.2))
axis(side = 2, at = seq(-0.4,0.4,0.2), lab = seq(-0.4,0.4,0.2), las = 2)

points(PCOA[,1], PCOA[,2],
       pch = metadat$shapes,
       cex = 1.8,
       xpd = NA,
       bg = add.alpha(metadat$colors, 0.7),
       col = add.alpha(metadat$colors, 0.9))
dev.off()

pdf(width = 8, height = 5, file = "pcoa_by_taxon_legend.pdf")
#Legend plot
plot(type = "n",
     PCOA[,1], PCOA[,2],
     ylim = c(0,10),
     axes = F,
     frame.plot = F,
     xlab = "",
     ylab = "",
     xlim = c(0,4)
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
#Things to do. 
#Figure out how to extract variance explained for each axis

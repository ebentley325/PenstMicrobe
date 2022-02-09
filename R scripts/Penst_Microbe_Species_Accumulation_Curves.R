#### Species accumulation Curves ####
#load packages
require(vegan)
#rownames(multi)<-multi$X
#multi<- multi[,-1]

#import the OTU table
otu16S<- read.csv("16S_otuTable_arranged_by_taxon.csv")
rownames(otu16S)<-otu16S$X
otu16S<- otu16S[,-1]
#otu16S[otu16S>0] <-1 #do this to change your data to prescence abscence data 
#totu16S<- t(otu16S)
par(mfrow = c(2,2))

#do the accumulation curve, random first
accumCurve<- specaccum(otu16S, method ="random", permutations = 1000)
#summary(accumCurve)
plot(accumCurve, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
#boxplot(accumCurve, col="yellow", add=TRUE, pch="+")


#do an accumulation curve, method = collector, which adds sites as they occur in the data. 
collectAccumCurve<-specaccum(otu16S, method ="collector")
plot(collectAccumCurve, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#use method = exact, which finds the expected (mean) species richness
exAccumCurve<- specaccum(otu16S, method ="exact")
plot(exAccumCurve, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#use method = rarefaction, finds the mean when accumulating individuals instead of sites

rareAccumeCurve<- specaccum(otu16S, method ="rarefaction")
plot(rareAccumeCurve, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")


mtext("Species Accumulation Curves 16S", side = 3, line = -21, outer = TRUE)




#### Species accumulation curves ITS ####
otuITS<- read.csv("ITS_otuTable_arranged_by_taxon.csv")
rownames(otuITS)<-otuITS$X
otuITS<- otuITS[,-1]
#otuITS[otuITS>0] <-1 #do this to change your data to prescence abscence data 
#totuITS<- t(otuITS)


#do the accumulation curve, random first
accumCurve<- specaccum(otuITS, method ="random", permutations = 1000)
#summary(accumCurve)
plot(accumCurve, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
#boxplot(accumCurve, col="yellow", add=TRUE, pch="+")


#do an accumulation curve, method = collector, which adds sites as they occur in the data. 
collectAccumCurve<-specaccum(otuITS, method ="collector")
plot(collectAccumCurve, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#use method = exact, which finds the expected (mean) species richness
exAccumCurve<- specaccum(otuITS, method ="exact")
plot(exAccumCurve, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#use method = rarefaction, finds the mean when accumulating individuals instead of sites

rareAccumCurve<- specaccum(otuITS, method ="rarefaction")
plot(rareAccumeCurve, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

mtext("Species Accumulation Curves ITS", side = 3, line = -31, outer = TRUE)



#### rarify the OTU table and re-run these things, control for sequencing depth ####
#16S
otuRare16S<-rrarefy(otu16S, 2000)
#totu16Srare<-t(otuRare16S)
accumCurveRare<- specaccum(otu16Srare, method ="random", permutations = 1000)
plot(accumCurveRare, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

collectAccumCurveRare<-specaccum(otu16Srare, method ="collector")
plot(collectAccumCurveRare, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

exAccumCurveRare<- specaccum(otu16Srare, method ="exact")
plot(exAccumCurveRare, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

rareAccumCurveRare<- specaccum(otu16Srare, method ="rarefaction")
plot(rareAccumCurveRare, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

mtext("Rarefied 16S Species Accumulation Curves ", side = 3, line = -21, outer = TRUE)

#ITS
otuRareITS<-rrarefy(otuITS, 50000)
#totuITSrare<-t(otuRareITS)
accumCurveRare<- specaccum(otuITSrare, method ="random", permutations = 1000)
plot(accumCurveRare, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

collectAccumCurveRare<-specaccum(otuITSrare, method ="collector")
plot(collectAccumCurveRare, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

exAccumCurveRare<- specaccum(otuITSrare, method ="exact")
plot(exAccumCurveRare, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

rareAccumCurveRare<- specaccum(otuITSrare, method ="rarefaction")
plot(rareAccumCurveRare, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

mtext("Rarefied ITS Species Accumulation Curves ", side = 3, line = -21, outer = TRUE)



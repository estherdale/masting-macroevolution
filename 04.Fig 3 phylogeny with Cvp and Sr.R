# plots phylogeny with rings of values of masting metrics and labels for plant families
# Written by Esther Dale, Feb 2021

library(phytools)
library(plotrix)

#### read in data and phylogeny ###
mast.df <- as.data.frame(readRDS("Masting data names checked 5 yrs.Rdata"))
rownames(mast.df) <- gsub(" ", "_", mast.df$Species)

phylo <- read.tree("Masting phylogeny 5yrs.nwk")

# look up higher taxonomic info #
mast.df$genus <- sapply(strsplit(mast.df$Species, split = " "), `[`, 1) #make genus column
tnrs.df <- read.csv("tnrs_genus.csv")
mast.df$family <- tnrs.df$Accepted_name_family[match(unlist(mast.df$genus), tnrs.df$Accepted_name)]
#identify clades
clade.list <-names( table(mast.df$family)[table(mast.df$family) >=8])
for(i in 1:length(clade.list)){
  rows.i <- which(mast.df$family==clade.list[i])
  clade.node <- getMRCA(phy=phylo, tip=(row.names(mast.df)[rows.i]))
  clade.row <- data.frame(cbind(clade=clade.list[i], mrca=clade.node))
  if(i==1){ clade.df <- clade.row
  }else{
  clade.df <- rbind(clade.df, clade.row)
}
}
#sort to same order as phylo
mast.df<- mast.df[match(unlist(phylo$tip.label), rownames(mast.df)),]

### generate colour values for each var for species ###
# all rescaled both were -1 to +1
# colours determined for each value
# for AR1 and synchrony white is 0 and blues are low masting and reds are extreme masting
rescale <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
ar1.rescale <- (mast.df$mean_AR1/2)+0.5
cv.rescale<- rescale(mast.df$mean_CV)
s.global.rescale <- (mast.df$mean_spat_global/2)+0.5
s.global.rescale[is.na(s.global.rescale)] <-0.5  #make NAs 0.5 so they will be white

color.ramp<- colorRamp(colors = c("#046C9A", "white", "#F2300F"), space="rgb")
color.ramp.inverse <- colorRamp(colors = c("#F2300F", "white", "#046C9A"), space="rgb")

ar1.cols <-rgb(color.ramp.inverse(ar1.rescale), maxColorValue = 255) 
cv.cols <-rgb(color.ramp(cv.rescale), maxColorValue = 255) 
s.global.cols <- rgb(color.ramp(s.global.rescale), maxColorValue = 255) 
s.global.cols[is.na(mast.df$mean_spat_global)] <- "black"

#make vector of width for rings on plot
set.width <- rep(20, nrow(mast.df))
names(set.width) <- rownames(mast.df)

#make shorter versions of long clade names
clade.names <- clade.list
clade.names[c(1,2,3,6,9)] <- gsub("aceae", "\\.", clade.names[c(1,2,3,6,9)])

# make ring arcs for the legend
legend.tips <- round(1.25*length(phylo$tip.label)/8, digits=0):round(2.75*length(phylo$tip.label)/8, digits = 0)
cv.legend <- rep(0, length(phylo$tip.label))
cv.legend[legend.tips] <- rgb(color.ramp(seq(from=1, to=0, length.out = length(legend.tips))), maxColorValue = 255) 
ar1.legend <- rep(0, length(phylo$tip.label))
ar1.legend[legend.tips] <- rgb(color.ramp(seq(from=1, to=0, length.out = length(legend.tips))), maxColorValue = 255) 
arrow.legend <-rep(0, length(phylo$tip.label))
arrow.legend[legend.tips[16:45]] <- "black"

pdf(file="Fig 2 phylo with masting metrics.pdf", height=5, width=5)
par(mar=c(0.5,0.5,0.5,0.5))
plot(phylo, type="fan", show.tip.label = T, x.lim=c(-100,100), label.offset =100 , tip.color = "white")
ring(set.width, phylo, col=ar1.cols, offset=5, cex=15, style="ring", lwd=5)
ring(set.width, phylo, col=cv.cols, offset=30, style="ring")
ring(set.width, phylo, col=s.global.cols, offset=55, cex=3, style="ring")
for(i in 1:length(clade.list)){
  arc.cladelabels(tree=NULL, text=clade.names[i], node=clade.df$mrca[i], mark.node=F, ln.offset = 1.275, lab.offset =1.35, lwd=2, cex=0.9)
}
#legend
ring(set.width, phylo, col=ar1.legend, offset=-200, style="ring")
ring(set.width, phylo, col=cv.legend, offset=-240, style="ring")
ring(set.width, phylo, col=ar1.legend, offset=-280, style="ring")
#arrow
ring(set.width/10, phylo, col=arrow.legend, offset=-155, style="ring")
arctext("weak <", center=c(0,0), radius=171, end =1.8, cex=0.75)
arctext("> strong", center=c(0,0), radius=171, start =1.33, cex=0.75)

#S scale
arctext("-1", center=c(0,0), radius=135, start =2.35, cex=0.6)
arctext("1", center=c(0,0), radius=135, start =0.93, cex=0.6)

#CVp scale
arctext("0.2", center=c(0,0), radius=95, start =2.5, cex=0.6)
arctext("2.7", center=c(0,0), radius=95, start =0.93, cex=0.6)

#ar1 scale
arctext("1", center=c(0,0), radius=55, start =2.4, cex=0.6)
arctext("-1", center=c(0,0), radius=55, start =0.93, cex=0.6)


arctext("S", center=c(0,0), radius=140, end =2.45, cex=0.7)
arctext("CVp", center=c(0,0), radius=100, end=2.6, cex=0.7)
arctext("AR1", center=c(0,0), radius=60, end=2.7, cex=0.7)

#scale bar
max.age <- max(nodeHeights(phylo))
axis(1, pos=-0.03*max.age,at=seq(50, 250, by=50),lwd=0.75, cex.axis=0.5, tck=-0.005, labels=F)
text(seq(50, 250, by=50), rep(-0.08*max.age, 5), seq(50, 250, by=50), cex=0.5)
text(x=0.5*max.age,y=-0.15*max.age,"Time (Ma)", cex=0.6)

dev.off()


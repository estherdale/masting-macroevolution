#makes Fig 4 of modal QuaSSE model for CVp and Synchrony

#read in data
duration=5
unit.type="full"
quasse.outputs <- readRDS(file=paste("Quasse_results_masting_", duration, "yrs_", unit.type, ".Rdata", sep=""))
models.i <- grep("modal", quasse.outputs$Model)


pdf(file="Quasse model plots Fig 4. modal CV S global.pdf", height=4, width=8)
par(mfrow=c(1,2), mar=c(4,4.1,1,0.5))
scale.factor <- 260 #scale diff betw freq and speciation axes

#CVp
curve(quasse.outputs$l.y0[models.i[2]] + ( quasse.outputs$l.y1[models.i[2]]-quasse.outputs$l.y0[models.i[2]])*exp(-(((x-quasse.outputs$l.xmid[models.i[2]])^2)/(2*quasse.outputs$l.s2[models.i[2]]))), from = 0,
      to=2.7, xlab = expression('CV'['p']), ylab="Speciation rate (Spp./Ma)", ylim=c(0,65), type="n", yaxt="n")
hist(quasse.cv, breaks=seq(from=0, by=0.195, length=15), add = T, border=NA, col="lightgrey")

axis(side = 2, at=seq(from=0, by=13, length=6), labels = seq(from=0, by=13/scale.factor, length=6))
curve(scale.factor*(quasse.outputs$l.y0[models.i[2]] + ( quasse.outputs$l.y1[models.i[2]]-quasse.outputs$l.y0[models.i[2]])*exp(-(((x-quasse.outputs$l.xmid[models.i[2]])^2)/(2*quasse.outputs$l.s2[models.i[2]])))), from = 0,
      to=2.7, lwd=1.5, add=T)
text(0.1, 0.24*scale.factor, "a)")

#global S
par(mar=c(4,0.5,1,4))
curve(quasse.outputs$l.y0[models.i[4]] + ( quasse.outputs$l.y1[models.i[4]]-quasse.outputs$l.y0[models.i[4]])*exp(-(((x-quasse.outputs$l.xmid[models.i[4]])^2)/(2*quasse.outputs$l.s2[models.i[4]]))), from = -1,
      to=1, xlab = expression(paste('S'['r'],' (Pearson correlation coefficient)'), sep=""), ylab="", ylim=c(0,65), type="n", yaxt="n")
hist(quasse.Sglobal, breaks=seq(from=-1, by=0.15, length=15), add = T, border=NA, col="lightgrey")
axis(side = 4, at=seq(from=0, by=10, length=7), labels = seq(from=0, by=10, length=7))
mtext(side=4, "Frequency", line=2.75)
curve(scale.factor*(quasse.outputs$l.y0[models.i[4]] + ( quasse.outputs$l.y1[models.i[4]]-quasse.outputs$l.y0[models.i[4]])*exp(-(((x-quasse.outputs$l.xmid[models.i[4]])^2)/(2*quasse.outputs$l.s2[models.i[4]])))), from = -1,
      to=1, lwd=1.5, add=T)
text(-0.92, 0.24*scale.factor, "b)")

dev.off()
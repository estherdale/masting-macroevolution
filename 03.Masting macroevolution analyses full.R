# Masting analyses of mastree data with different duration thresholds and unit type combinations
# 1. Phylogenetic signal and correlation
# 2. Trait-dependent evolution
# 3. Compares macroevolutionary rates of MASTREE+ data to other species in the Smith and Brown GBOTB phylogeny
# 4. Cluster analysis based on masting matrics
# 5. Compare synchrony, latitude, and macro rates between clusters
# Written by Esther Dale Feb 2021

#required packages and functions

library(ape)
library(BAMMtools)
library(caper)
library(cluster)
library(DDD)
library(diversitree)
library(dplyr)
library(geiger)
library(phytools)
library(picante)
library(phylolm)
library(secsse)
library(tidyr)

source("AIC function.R")
source("DRmetric.R")
source("essim function.R")
source("Tip ages function.R")
source("Plot text rounding functions.R")


###############################################
############### prep files#####################
###############################################
masting.datasets <-list.files(pattern="Masting data names checked .") #masting data for different durations/unit types
masting.phylos <-list.files(pattern="Masting phylogeny .")

#load BAMM object for Smith and Brown phylogeny
Smith.Brown.BAMM <- readRDS("SmithBrown_BAMM_vascular.RDS")
#BAMM object from Igea & Tanentzap 2020 https://github.com/javierigea/LDGplants_rates

#df with original tip label info so can match tips in SmithBrown.BAMM
tip.labs.checked <- readRDS("TPL output GBOTB names masting genera.Rdata")
  #tip.labs.checked$Taxon are the original names

# seed size evolution BAMM
seed.BAMM.full <- readRDS("Trait_edata_50shifts_1000samples.Rsave")
#BAMM object from Igea et al. 2017 https://github.com/javierigea/seed_size
seed.names<- readRDS("Seed names checked.Rdata")
  #seed.names$Taxon are the original names from the bamm object

# Smith and Brown GBOTB phylogeny
Smith.brown.tree <- read.tree("GBOTB.tre")
#phylogeny from Smith & Brown 2018 https://doi.org/10.1002/ajb2.1019 
#tree file available here https://github.com/FePhyFoFum/big_seed_plant_trees

# DR values for Smith and Brown phylogeny
  #this takes 5+ hours
smith.brownDR <- DR_statistic(Smith.brown.tree)

# prep age, lambda, beta, and dr rates for tips
tip.evolution.df <- calculate_tip_ages(Smith.brown.tree)
bamm.rates <- getTipRates(Smith.Brown.BAMM)
tip.evolution.df$lambda <- bamm.rates$lambda.avg[match(unlist(tip.evolution.df$tip), names(bamm.rates$lambda.avg))]
seed.rates <- getTipRates(seed.BAMM.full)
tip.evolution.df$beta <- seed.rates$beta.avg[match(unlist(tip.evolution.df$tip), names(seed.rates$beta.avg))]
tip.evolution.df$dr <- smith.brownDR[match(unlist(tip.evolution.df$tip), names(smith.brownDR))]
tip.evolution.df$tip.age <- as.numeric(as.character(tip.evolution.df$tip.age))
rownames(tip.evolution.df) <- tip.evolution.df$tip
tip.evolution.df$tip <- as.character(tip.evolution.df$tip)

###################################################
###################################################
for(j in 1:length(masting.datasets)){

mast.df <- as.data.frame(readRDS(masting.datasets[j]))
rownames(mast.df) <- gsub(" ", "_", mast.df$Species)

duration <- as.numeric(substr(masting.datasets[j], 28, 28)) #number of years minimum duration threshold
unit.type <- rep(c("area", "count", "individual", "mass", "full"),4)[j]

phylo <- read.tree(masting.phylos[duration-2])
phylo <- keep.tip(phylo, tip = unique(gsub(" ", "_", mast.df$Species)))
phylo$node.label <- NA

################################################
########### 1. Phylogenetic signal #############
################################################

 #vars to have phylogenetic signal tested
 vars <- colnames(mast.df)[c(6,15:16)]

 #sets up empty dataframe for phylogenetic signal results
   # var=variable, K.stat= Blomberg's K statistic, K.sd=standard deviation of K, K.p=p value testing significance of k
 K.results <- data.frame("Var"=character(), "Method"=character(),"Stat.without.se"=numeric(), "Stat.with.se"=numeric(), "sigma.sq"=numeric(),
                         "p.without.se"=numeric(), "p.with.se"=numeric(), "n.without.se"=numeric(), "n.with.se"=numeric(), stringsAsFactors = F)

 # loops through variables
 for(i in 1:length(vars)){
   rows.i <- (2*i-1):(2*i)
   col.sd <- grep(gsub("mean", "sd", vars[i]), colnames(mast.df))[1] #std dev column for var i

   #full data
   mast.i <- mast.df[,c(1, grep(vars[i], colnames(mast.df))[1])]
   mast.i$Species <- gsub(" ", "_", mast.i$Species)
   phylo.i <- keep.tip(phylo, tip=mast.i$Species)
   mast.cd <- comparative.data(phylo.i, mast.i, names.col=Species)

   #for data with std dev - want only species with std dev value >0 and std dev not more than 5x higher than mean
   rows.sd <- complete.cases(mast.df[,col.sd])&mast.df[,col.sd]>0 & (abs(mast.df[,col.sd]/mast.i[,2]))<5
   mast.i.se <- mast.df[rows.sd,c(1, grep(vars[i], colnames(mast.df))[1], col.sd, 8)]
   mast.i.se$Species <-  gsub(" ", "_", mast.i.se$Species)
   #make column with std error
   mast.i.se$se <- mast.i.se[,3]/sqrt(mast.i.se$n_timeseries)
   phylo.se <-  keep.tip(phylo, tip=mast.i.se$Species)
   phylo.se$node.label <- NA
   mast.se.cd <- comparative.data(phylo.se, mast.i.se, names.col = Species)

   #k with and without std err
   K.i.se <- phylosig(mast.se.cd$phy, mast.se.cd$data[,1],se=mast.se.cd$data$se, method="K", test = T)
   K.i <- phylosig(mast.se.cd$phy, mast.se.cd$data[,1], method="K", test = T)

   # lambda with and without standard error
       #if model throws an error then try different starting values of lambda, and if that doesn't work, different fixed values of lambda and picks best option
   yy <- mast.se.cd$data[,1]
   y2 <- mast.cd$data[,1]
   ee <- 1/mast.se.cd$data$se

   lamb.se <- tryCatch({ gls(yy ~ 1, correlation=corPagel(0.5, mast.se.cd$phy,fixed=F), method="ML", weights=varFixed(~ee)) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
   if(is.null(lamb.se)){
     lamb.se <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(yy ~ 1, correlation=corPagel(y, mast.se.cd$phy,fixed=F), method="ML", weights=varFixed(~ee))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
     if(all(sapply(lamb.se, is.null) == T)){
       lamb.se <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(yy ~ 1, correlation=corPagel(y, mast.se.cd$phy,fixed=T), method="ML", weights=varFixed(~ee)) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
     }
     lamb.se <- lamb.se[which(sapply(lamb.se, is.null) == F)][[which.min(sapply(lamb.se[which(sapply(lamb.se, is.null) == F)],AIC))]]
   }
   try(if(summary(lamb.se)$modelStruct[1]<0){
    lamb.se <- try(gls(yy ~ 1, correlation=corPagel(0,mast.se.cd$phy,fixed=T), method="ML"))
   })
   lamb.null.se <- gls(yy~1, method="ML", weights=varFixed(~ee))
   lamb.test.se <- try(anova(lamb.null.se, lamb.se))
   lamb.se.p <- try(lamb.test.se$`p-value`[2])
   if(is.null(lamb.se.p)){lamb.se.p <- 1}

  lamb.i <- tryCatch({ gls(y2 ~ 1, correlation=corPagel(0.5, mast.cd$phy,fixed=F), method="ML" ) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  if(is.null(lamb.i)){
    lamb.i <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(y2 ~ 1, correlation=corPagel(y, mast.cd$phy,fixed=F), method="ML")},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
    if(all(sapply(lamb.i, is.null) == T)){
      lamb.i <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(y2 ~ 1, correlation=corPagel(y, mast.cd$phy,fixed=T), method="ML") }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
    }
    lamb.i <- lamb.i[which(sapply(lamb.i, is.null) == F)][[which.min(sapply(lamb.i[which(sapply(lamb.i, is.null) == F)],AIC))]]
  }

   try(if(summary(lamb.i)$modelStruct[1]<0){
       lamb.i <-try( gls(y2 ~ 1, correlation=corPagel(0,mast.cd$phy,fixed=T), method="ML"))
       })
  lamb.null <- gls(y2~1, method="ML")
   lamb.test <- try(anova(lamb.null, lamb.i))
   lamb.p <- try(lamb.test$`p-value`[2])
   if(is.null(lamb.p)){lamb.p <- 1}

#fill in dataframe
  K.results[rows.i[1],] <- c(vars[i],"K", K.i$K, K.i.se$K, K.i.se$sig2, K.i$P, K.i.se$P, length(mast.cd$phy$tip.label), length(mast.se.cd$phy$tip.label))
   K.results[rows.i[2],] <- c(vars[i],"Lambda", NA,NA, NA, NA, NA, length(mast.cd$phy$tip.label), length(mast.se.cd$phy$tip.label))
  try( K.results[rows.i[2],c( "Stat.without.se","p.without.se" )] <- c(summary(lamb.i)$modelStruct[1],lamb.p))
   try( K.results[rows.i[2],c( "Stat.with.se","p.with.se" )] <- c(summary(lamb.se)$modelStruct[1],lamb.se.p))

   try(rm(lamb.se))
 try(rm(lamb.i))
   }

 K.results$dataset <- rep(duration, nrow(K.results))
 K.results$unit.type <- rep(unit.type, nrow(K.results))

 saveRDS(K.results, file=paste("Phylogenetic_signal_results_", duration, "yrs_", unit.type, ".Rdata", sep=""))

 ## phylogenetic correlation ##
 corr.results.df <- data.frame("vars"=character(), "R-squared"=numeric(),"r"=numeric(), "coefficient"=numeric(), "p-value"=numeric(), "n"=numeric())
 mast.cors <- list(c("mean_CV", "mean_spat_global"), c("mean_CV", "mean_spat_100"), c("mean_spat_global", "mean_spat_100"))
 for(i in 1:length(mast.cors)){
   #columns of vars to be correlated
   cols.i <-c( grep(mast.cors[[i]][1], colnames(mast.df)), grep(mast.cors[[i]][2], colnames(mast.df)))

   #set up comparative data object
   cor.df <- mast.df[,c(1, cols.i)]
   cor.df <- cor.df[complete.cases(cor.df),]
   cor.df$Species <- gsub(" ", "_", cor.df$Species)
   colnames(cor.df)[2:3] <- c("x", "y")
   phylo.i <- keep.tip(phylo, cor.df$Species)
   phylo.i$node.label <- NA
   cor.cd <- comparative.data(phylo.i, cor.df, names.col = 'Species')
   #fit pgls
   try(cor.pgls.i <- pgls(y~x, cor.cd, lambda="ML"))
   try(cor.summary.i <- summary(cor.pgls.i))
   #calculate r
   try(r.i <- sqrt(cor.summary.i$r.squared))
   try(if(cor.summary.i$coefficients[2,1]<0){r.i <- r.i*-1})

   #put results in table
   corr.results.df[i,] <- c(paste(mast.cors[[i]][1], mast.cors[[i]][2], sep=" & "), rep(NA, 4), nrow(cor.df))
   try(corr.results.df[i,2:5] <- c(cor.summary.i$r.squared, r.i, cor.summary.i$coefficients[2,1], cor.summary.i$coefficients[2,4]))

  rm(list=c('cor.pgls.i', 'cor.summary.i'))
}

corr.results.df$dataset <- rep(duration, nrow(corr.results.df))
corr.results.df$unit.type <- rep(unit.type, nrow(corr.results.df))

 saveRDS(corr.results.df, file=paste("Phylogenetic_correlation_results_", duration, "yrs_", unit.type, ".Rdata", sep=""))

######################################################
########### 2. Trait dependent evolution #############
######################################################

## Analyses to test relationship between masting traits (cv, ar1 and synchrony) and speciation and seed size evolution
# uses STRAPP, ESSIM, and QuaSSE approaches

#### prep BAMM objects#####
## diversification ##
#get original phylo tip labels and add "_" so can match to BAMM object
mast.df$tips_original <- tip.labs.checked$Taxon[match(unlist(mast.df$Species), paste(tip.labs.checked$New.Genus, tip.labs.checked$New.Species, sep=" "))]
mast.df$tips_original <- gsub(" ", "_", mast.df$tips_original)
row.names(mast.df) <- mast.df$tips_original

# subset bamm object to species in masting
mast.BAMM <- subtreeBAMM(Smith.Brown.BAMM, tips = mast.df$tips_original, node = NULL)
mast.bamm.df <- mast.df[mast.df$tips_original %in% mast.BAMM$tip.label,]

# drop NAs from tip states
mast.BAMM$tipLambda <- 
              mast.BAMM$tipLambda[which(sapply(mast.BAMM$tipStates,
                                                        function(x){any(is.na(x))==F}))]
mast.BAMM$tipStates <-  
              mast.BAMM$tipStates[which(sapply(mast.BAMM$tipStates,
                                   function(x){any(is.na(x))==F}))]

## seed size evolution ##
# get the masting data to match the bamm object
seed.df <- mast.df #copy the masting data
seed.df$tips_original <- seed.names$Taxon[match(unlist(mast.df$Species), seed.names$binomial)]
seed.df$tips_original <- gsub(" ", "_", seed.df$tips_original)
seed.df <- seed.df[complete.cases(seed.df$tips_original),]
rownames(seed.df) <- seed.df$tips_original

# subset bamm object to species in masting
seed.BAMM <- subtreeBAMM(seed.BAMM.full, tips = seed.df$tips_original, node = NULL)
seed.df <- seed.df[seed.df$tips_original %in% seed.BAMM$tip.label,]

#### prep for ES-sim analysis #####
#subset to masting species
mast.DR <- smith.brownDR[names(smith.brownDR) %in% mast.df$tips_original] #subset to masting species

#### prep for QuaSSE ####
#set up data vectors for traits without missing values
quasse.ar1 <- mast.df$mean_AR1[complete.cases(mast.df$mean_AR1)]
quasse.cv <- mast.df$mean_CV[complete.cases(mast.df$mean_CV)]
quasse.S100 <- mast.df$mean_spat_100[complete.cases(mast.df$mean_spat_100)]
quasse.Sglobal <- mast.df$mean_spat_global[complete.cases(mast.df$mean_spat_global)]
#add names
names(quasse.ar1) <- gsub(" ", "_", mast.df$Species[complete.cases(mast.df$mean_AR1)])
names(quasse.cv) <- gsub(" ", "_", mast.df$Species[complete.cases(mast.df$mean_CV)])
names(quasse.S100) <- gsub(" ", "_", mast.df$Species[complete.cases(mast.df$mean_spat_100)])
names(quasse.Sglobal) <- gsub(" ", "_", mast.df$Species[complete.cases(mast.df$mean_spat_global)])
#make into list
quasse.states <- list(quasse.cv, quasse.S100, quasse.Sglobal)

#set up standard error vectors for quasse mean estimates
quasse.se.ar1 <- mast.df$sd_AR1[complete.cases(mast.df$mean_AR1)]/sqrt(mast.df$n_timeseries[complete.cases(mast.df$mean_AR1)])
quasse.se.cv <- mast.df$sd_CV[complete.cases(mast.df$mean_CV)]/sqrt(mast.df$n_timeseries[complete.cases(mast.df$mean_CV)])
quasse.se.S100 <- mast.df$sd_spat_100[complete.cases(mast.df$mean_spat_100)]/sqrt(mast.df$N_spat_100[complete.cases(mast.df$mean_spat_100)])
quasse.se.Sglobal <- mast.df$sd_spat_global[complete.cases(mast.df$mean_spat_global)]/sqrt(mast.df$N_spat_global[complete.cases(mast.df$mean_spat_global)])

#fill in 0s/NAs with mean value
quasse.se.ar1[quasse.se.ar1==0] <- mean(quasse.se.ar1[quasse.se.ar1>0])
quasse.se.cv[quasse.se.cv==0] <- mean(quasse.se.cv[quasse.se.cv>0])
quasse.se.S100[is.na(quasse.se.S100)|quasse.se.S100==0] <- mean(quasse.se.S100, na.rm = T)
quasse.se.Sglobal[is.na(quasse.se.Sglobal)|quasse.se.Sglobal==0] <- mean(quasse.se.Sglobal, na.rm = T)

#add names
names(quasse.se.ar1) <- gsub(" ", "_", mast.df$Species[complete.cases(mast.df$mean_AR1)])
names(quasse.se.cv) <- gsub(" ", "_", mast.df$Species[complete.cases(mast.df$mean_CV)])
names(quasse.se.S100) <- gsub(" ", "_", mast.df$Species[complete.cases(mast.df$mean_spat_100)])
names(quasse.se.Sglobal) <- gsub(" ", "_", mast.df$Species[complete.cases(mast.df$mean_spat_global)])
#make into list
quasse.se.states <- list(quasse.se.cv, quasse.se.S100, quasse.se.Sglobal)

#prune phylogenies and force ultrametric
phylo.quasse <- keep.tip(phylo, names(quasse.ar1))
phylo.quasse <- force.ultrametric(phylo.quasse)
phylo.S100 <- keep.tip(phylo, names(quasse.S100))
phylo.S100 <- force.ultrametric(phylo.S100 )
phylo.Sglobal <- keep.tip(phylo, names(quasse.Sglobal))
phylo.Sglobal <- force.ultrametric(phylo.Sglobal)

#list of the phylogenies
quasse.phylos <- list(phylo.quasse, phylo.S100, phylo.Sglobal)

##### testing associations between speciation/seed evolution and masting traits ####
#set up dataframes
trait.divn.results <- data.frame("Variable"=character(), "Evolution.type"=character(),"Method"=character(), "rho"=numeric(), "p"=numeric(), "n"=numeric(), stringsAsFactors = F)
quasse.outputs <- data.frame("Variable"=character(), "Model"=character(), "lambda"=numeric(), "l.m"=numeric(),
                             "l.y0"=numeric(), "l.y1"=numeric(), "l.xmid"=numeric(), "l.r"=numeric(), "l.s2"=numeric(),
                             "mu"=numeric(), "diffusion"=numeric(), "n"=numeric(), "lnLik"=numeric(), "AIC"=numeric(),
                             "AIC change"=numeric())

#loop through traits
for(i in 1:length(vars)){
  col.i <- grep(vars[i], colnames(mast.df))[1]
  rows.i <- ((3*i)-2):(3*i) #rows for trait.divn.results df
  rows.iq <- ((4*i)-3):(4*i) #rows for quasse.outputs df
  
  # speciation STRAPP
  mast.i <- mast.bamm.df[,col.i]
  names(mast.i) <- row.names(mast.bamm.df)
  strapp.mast.i <- traitDependentBAMM(mast.BAMM, mast.i, 1000, rate = "speciation", method = "spearman", logrates = F ,two.tailed = TRUE, traitorder = "p", nthreads = 1)
  trait.divn.results[rows.i[1],] <- c(vars[i], "Speciation", "STRAPP", strapp.mast.i$estimate, strapp.mast.i$p.value, length(mast.i[complete.cases(mast.i)]))

  # seed evolution STRAPP
  seed.i <- seed.df[,col.i]
  names(seed.i) <- row.names(seed.df)
 strapp.seed.i <- traitDependentBAMM(seed.BAMM, seed.i, 1000, rate = "trait", method = "spearman", logrates = F ,two.tailed = TRUE, traitorder = "p", nthreads = 1)
  trait.divn.results[rows.i[2],] <- c(vars[i], "Trait", "STRAPP", strapp.seed.i$estimate, strapp.seed.i$p.value, length(seed.i[complete.cases(seed.i)]))

  # speciation ES-sim
  mast.ib <- mast.df[,col.i]
  names(mast.ib) <- gsub(" ", "_", mast.df$tips_original) #mast.ib has original names
  mast.ib <- mast.ib[complete.cases(mast.ib)]
  mast.tree.i <- keep.tip(Smith.brown.tree, names(mast.ib)) #Smith brown tree has unchecked names (unlike "phylo)
  mast.DR.i <- mast.DR[names(mast.DR) %in% names(mast.ib)]
  essim.mast.i<- essim(mast.tree.i, mast.ib, nsim=1000, es=mast.DR.i)
  trait.divn.results[rows.i[3],] <- c(vars[i], "Speciation","ES-sim", essim.mast.i[1], essim.mast.i[2], length(mast.ib[complete.cases(mast.ib)]))

  # diversification Quasse
  q.states.i <- quasse.states[[i]]
  q.se.i <- quasse.se.states[[i]]
  phylo.q.i <- quasse.phylos[[i]]
  
  #set up speciation functions (constant extinction for all)
  make.models <- function(lambda, mu) make.quasse(phylo.q.i, q.states.i, q.se.i, lambda, mu)
  nodrift <- function(f) constrain(f, drift ~ 0)
  p <- starting.point.quasse(phylo.q.i, q.states.i)
  xr <- range(q.states.i) + c(-1,1) * 20 * p["diffusion"]
  linear.x <- make.linear.x(xr[1], xr[2])
  
  quasse.const <- make.models(constant.x, constant.x) #constant speciation
  quasse.lin <- make.models(linear.x, constant.x) #linear function (within a range) speciation
  quasse.sig <- make.models(sigmoid.x, constant.x) #sigmoid speciation
  quasse.modal <- make.models(noroptimal.x, constant.x) #sigmoid speciation
  
  control <- list(parscale=.1, reltol=0.001)
  #fit constant model
  try(  mle.const <- find.mle(nodrift(quasse.const), p, lower=0, control=control,verbose=0))
  #set up starting vlaues for other models
  p.c <- mle.const$par
  p.l <- c(p.c[1], l.m=0, p.c[2:3])
  p.s <- p.m <- c(p.c[1], p.c[1], mean(xr),1,p.c[2:3])
  names(p.s) <- argnames(nodrift(quasse.sig))
  names(p.m) <- argnames(nodrift(quasse.modal))
  #fit linear, sigmoidal, and modal models
  try( mle.lin <- find.mle(nodrift(quasse.lin), p.l, control=control, verbose=0) )
  try( mle.sig <- find.mle(nodrift(quasse.sig), p.s, control=control, verbose=0))
  try( mle.modal <- find.mle(nodrift(quasse.modal), p.m, control=control, verbose=0))
  #compare models
  model.comp <- anova(mle.const, linear=mle.lin, sigmoidal=mle.sig, modal=mle.modal)
 
  # fill in results
  quasse.outputs[rows.iq[1],] <- c(vars[i], "constant", mle.const$par[1], rep(NA, 6), mle.const$par[2:3], length(q.states.i), model.comp[1,2:3], model.comp[1,3]-min(model.comp[,3]))
  quasse.outputs[rows.iq[2],] <- c(vars[i], "linear", mle.lin$par[1:2], rep(NA, 5), mle.lin$par[3:4], length(q.states.i), model.comp[2,2:3], model.comp[2,3]-min(model.comp[,3]), NA)
  quasse.outputs[rows.iq[3],] <- c(vars[i], "sigmoidal", rep(NA,2), mle.sig$par[1:4], NA, mle.sig$par[5:6], length(q.states.i), model.comp[3,2:3], model.comp[3,3]-min(model.comp[,3]))
  quasse.outputs[rows.iq[4],] <- c(vars[i], "modal", rep(NA,2), mle.modal$par[1:3], NA, mle.modal$par[4:6], length(q.states.i), model.comp[4,2:3], model.comp[4,3]-min(model.comp[,3]))
  
}

#add dataset info to dataframes
trait.divn.results$dataset <- rep(duration, nrow(trait.divn.results))
trait.divn.results$unit.type <- rep(unit.type, nrow(trait.divn.results))
quasse.outputs$dataset <- rep(duration, nrow(quasse.outputs))
quasse.outputs$unit.type <- rep(unit.type, nrow(quasse.outputs))

saveRDS(trait.divn.results, file=paste("Trait_dependent_diversification_results_masting_", duration, "yrs_", unit.type, ".Rdata", sep=""))
saveRDS(quasse.outputs, file=paste("Quasse_results_masting_", duration, "yrs_", unit.type, ".Rdata", sep=""))

########################################################################
########### 3. Compare masting data to non-mastree species #############
########################################################################
#make column for masting species in tip evolution dataframe
 tip.evolution.df$masting <- rep(0, nrow(tip.evolution.df))
 tip.evolution.df$masting[tip.evolution.df$tip %in% mast.df$tips_original] <- 1
 tip.evolution.df$masting <- as.factor(tip.evolution.df$masting)

# drop NAs from tip states
Smith.Brown.BAMM$tipLambda <-
  Smith.Brown.BAMM$tipLambda[which(sapply(Smith.Brown.BAMM$tipStates,
                                   function(x){any(is.na(x))==F}))]
Smith.Brown.BAMM$tipStates <-
  Smith.Brown.BAMM$tipStates[which(sapply(Smith.Brown.BAMM$tipStates,
                                   function(x){any(is.na(x))==F}))]

# strapp masting and lambda
# speciation STRAPP
lambda.j <- tip.evolution.df$masting
names(lambda.j) <- tip.evolution.df$tip
lambda.j <- lambda.j[complete.cases(tip.evolution.df$lambda)]

  tip.evolution.df[tip.evolution.df$tip %in% Smith.Brown.BAMM$tip.label,]
Smith.Brown.BAMM.j <- subtreeBAMM(Smith.Brown.BAMM, tips = lambda.j, node = NULL)

strapp.lambda.j <- traitDependentBAMM(Smith.Brown.BAMM.j, lambda.j, 1,
                                      rate = "speciation", method = "mann-whitney", logrates = F ,two.tailed =
                                        TRUE, traitorder = "1")

# seed evolution STRAPP
# subset bamm object to species in masting
seed.j <- tip.evolution.df$masting
names(seed.j) <- tip.evolution.df$tip
seed.j <- seed.j[names(seed.j) %in% seed.BAMM.full$tip.label]
seed.BAMM.j <- subtreeBAMM(seed.BAMM.full, tips = names(seed.j), node = NULL)

seed.j <- seed.j[names(seed.j)%in% seed.BAMM.j$tip.label]
strapp.seed.j <- traitDependentBAMM(seed.BAMM.j, seed.j, 1000, rate = "trait", method = "mann-whitney", logrates = F ,two.tailed = TRUE, traitorder = c("0","1"), nthreads = 1)

## phylo glms for age and dr
#rescale variables and log transform for using in phyloglm
 tip.evolution.df$age.scaled <- scale(log(tip.evolution.df$tip.age))
 tip.evolution.df$dr.scaled <- scale(log(tip.evolution.df$dr))

dr.glm1 <- phyloglm(masting~dr.scaled, phy=Smith.brown.tree, data=tip.evolution.df, method="logistic_MPLE", start.alpha = 0.0001, log.alpha.bound =6, btol = 30)

#put results in dataframe
bias.df <- data.frame("Variable"=character(), "method"=character(), "coefficient"=numeric(), "p-value"=numeric(), "median.mastree"=numeric(), "median.non-mastree"=numeric(), "N"=numeric())

bias.df[1,] <- c("lambda", "STRAPP", NA, strapp.lambda.j$p.value, unlist(strapp.lambda.j$estimate), length(lambda.j))
bias.df[2,] <- c("seed mass change", "STRAPP", NA, strapp.seed.j$p.value, unlist(strapp.seed.j$estimate), length(seed.j))
bias.df[3,] <- c("dr", "phyloglm", summary(dr.glm1)$coefficients[2,1], summary(dr.glm1)$coefficients[2,4], NA, NA, nrow(tip.evolution.df))

bias.df$dataset <- rep(duration, nrow(bias.df))
bias.df$unit.type <- rep(unit.type, nrow(bias.df))

saveRDS(bias.df, file=paste("Macroevolutionary_bias_mastree_data_", duration, "yrs_", unit.type, ".Rdata", sep=""))

### plots ###
masting.colours <-c("#046C9A", "#F2300F")

pdf(file=paste("Fig S4 boxplots of mastree and non-mastree macroevolutionary rates ", duration,"yrs ", unit.type, ".pdf", sep=""), height = 6, width = 6)
par(mfrow=c(2,2), mar=c(4,4.1,1,1.5))

# lambda
boxplot(lambda~masting, data=tip.evolution.df, main="", ylab="Rate of speciation (spp./Ma)", log="y", xaxt="n", yaxt="n", col=masting.colours, xlab="")
title(xlab="Species", line=2.25)
axis(1,at=c(1,2), labels=c("non-MASTREE+", "MASTREE+"), cex.axis=0.8)
axis(2, at=c(0.005, 0.05, 0.5, 5), labels=c(0.005, 0.05, 0.5, 5))
text(0.45, exp(0.95*max(log(tip.evolution.df$lambda), na.rm = T)), "a)", pos=4)
text(1.2,exp(0.9*max(log(tip.evolution.df$lambda), na.rm = T)), "Median rates:", pos=4, cex=0.75)
text(1.3,exp(0.75*max(log(tip.evolution.df$lambda), na.rm = T)), paste("non-MASTREE+ =", round(as.numeric(bias.df$median.non-mastree[2]), digits = 2)), pos=4, cex=0.75)
text(1.3,exp(0.60*max(log(tip.evolution.df$lambda), na.rm = T)), paste("MASTREE+ =", round(as.numeric(bias.df$median.mastree[2]), digits = 2)), pos=4, cex=0.75)
text(1.2,exp(0.45*max(log(tip.evolution.df$lambda), na.rm = T)), paste("p-value =", round(as.numeric(bias.df$p.value[2]), digits = 2)), pos=4, cex=0.75)

# seed mass change
boxplot(beta~masting, data=tip.evolution.df, main="", ylab="Rate of seed mass change", log="y",
        xaxt="n", yaxt="n", col=masting.colours, xlab="")
title(xlab="Species", line=2.25)
axis(1,at=c(1,2), labels=c("non-MASTREE+", "MASTREE+"), cex.axis=0.8)
axis(2, at=c(0.005, 0.05, 0.5, 5), labels=c(0.005, 0.05, 0.5, 5))
text(0.45, exp(0.95*max(log(tip.evolution.df$beta), na.rm = T)), "b)", pos=4)
text(1.2,exp(0.9*max(log(tip.evolution.df$beta), na.rm = T)), "Median rates:", pos=4, cex=0.75)
text(1.3,exp(0.73*max(log(tip.evolution.df$beta), na.rm = T)), paste("non-MASTREE+ =", round(as.numeric(bias.df$median.non-mastree[3]), digits = 2)), pos=4, cex=0.75)
text(1.3,exp(0.56*max(log(tip.evolution.df$beta), na.rm = T)), paste("MASTREE+ =", round(as.numeric(bias.df$median.mastree[3]), digits = 2)), pos=4, cex=0.75)
text(1.2,exp(0.39*max(log(tip.evolution.df$beta), na.rm = T)), paste("p-value =", round(as.numeric(bias.df$p.value[3]), digits = 2)), pos=4, cex=0.75)

# DR
boxplot(dr~masting, data=tip.evolution.df, main="", ylab="Diversification rate metric (DR, spp./Ma)", log="y",
        xaxt="n", yaxt="n",  col=masting.colours, xlab="")
title(xlab="Species", line=2.25)
axis(1,at=c(1,2), labels=c("non-MASTREE+", "MASTREE+"), cex.axis=0.8)
axis(2, at=c(0.01, 0.1, 1, 10, 100), labels=c(0.01, 0.1, 1, 10, 100))
text(0.45, exp(0.95*max(log(tip.evolution.df$dr), na.rm = T)), "c)", pos=4)
text(1.5,exp(0.9*max(log(tip.evolution.df$dr), na.rm = T)), paste("Coefficient =", round(as.numeric(bias.df$coefficient[4]), digits = 2)), pos=4, cex=0.75)
text(1.5,exp(0.8*max(log(tip.evolution.df$dr), na.rm = T)), "p-value < 0.001", pos=4, cex=0.75)

dev.off()

######################################################
############### 4. Cluster analysis #####################
######################################################
#generate dissimilarity matrix using Gower's distance
mast.clust <- daisy(mast.df[,c(4,6)], metric="gower")

#determine within-cluster sum of squares
WCSS <- c()
for(i in 2:12){
  km <- kmeans(mast.clust, i)
  WCSS <- c(WCSS, mean(km$withinss))
}

#plot to identify elbow
plot.name <- paste("Elbow plot masting clusters ", duration, "yrs ", unit.type, ".pdf", sep="")
pdf(file=plot.name)
plot(2:12, WCSS, pch=19, xlab="Number of clusters", ylab="Mean sum of squares")
lines(2:12, WCSS)
dev.off()

# will use 5 clusters going forward
no.clusters <- 5
mast.clust <- kmeans(mast.clust, no.clusters)

#add cluster to masting data
mast.df$cluster <- as.factor(mast.clust$cluster)

#add in clusters
tip.evolution.df$cluster <- rep(0, nrow(tip.evolution.df))
for(i in 1:no.clusters){
  tip.evolution.df$cluster[tip.evolution.df$tip %in% mast.df$tips_original[mast.df$cluster==i]] <- i
}
tip.evolution.df$cluster <- as.factor(tip.evolution.df$cluster)

##################################################################
################# 5. Comparing clusters ##########################
##################################################################
## set up SsecSSE
clusters.j <- data.frame(gsub(" ", "_", mast.df$Species), mast.df$cluster)
phylo.clusters <- force.ultrametric(phylo)
colnames(clusters.j) <- c("Species", "Cluster")
clusters.j <- sortingtraits(clusters.j, phylo.clusters)
c.states <- 5 #number of concealed states (should match number of examined states)

#set up pars list
idparslist.cr <- id_paramPos(clusters.j, num_concealed_states = c.states)
idparslist.etd <- idparslist.cr
idparslist.ctd <- idparslist.cr

#set starting number for first transition rate parameter 
# =number of speciation params + number of extinction parameters +1)
k.values <- c(3,7,7)
count.cr <- k.values[1]
count.etd <- k.values[2]
count.ctd <- k.values[3]

#fill in diagonals as NAs and dual transitions (both visible state and concealed state transition) to 0
# all transitions set to the same rate
for(i in 0:((c.states^2)-1)){
  for(h in 0:((c.states^2)-1)){
    if(i==h){  
      idparslist.cr[[3]][(i+1),(h+1)] <- NA 
      idparslist.etd[[3]][(i+1),(h+1)] <- NA 
      idparslist.ctd[[3]][(i+1),(h+1)] <- NA 
    }else{
      if(i%/%c.states==h%/%c.states | i%%c.states==h%%c.states){
        idparslist.cr[[3]][(i+1),(h+1)] <- count.cr
        idparslist.etd[[3]][(i+1),(h+1)] <- count.etd
        idparslist.ctd[[3]][(i+1),(h+1)] <- count.ctd
        #  count <- count+1 #have this if want to have different transition rates
      }else{
        idparslist.cr[[3]][(i+1),(h+1)] <- 0
        idparslist.etd[[3]][(i+1),(h+1)] <- 0
        idparslist.ctd[[3]][(i+1),(h+1)] <- 0
      }
    }
    
  }#h
} #i

#make all extinction rates the same (=number of speciation parameters+1)
idparslist.cr[[2]][] <- 2
idparslist.etd[[2]][] <- 6
idparslist.ctd[[2]][] <- 6

#set up different speciation options
#CR model all speciation rates the same
idparslist.cr[[1]][1:25] <- 1
#ETD model all speciation rates the same within examined states (numbered)
idparslist.etd[[1]][c(1,6,11,16,21)] <- 1
idparslist.etd[[1]][c(2,7,12,17,22)] <- 2
idparslist.etd[[1]][c(3,8,13,18,23)] <- 3
idparslist.etd[[1]][c(4,9,14,19,24)] <- 4
idparslist.etd[[1]][c(5,10,15,20,25)] <- 5
#CTD model all speciation rates the same within concealed states (lettered)
idparslist.ctd[[1]][1:5] <- 1
idparslist.ctd[[1]][6:10] <- 2
idparslist.ctd[[1]][11:15] <- 3
idparslist.ctd[[1]][16:20] <- 4
idparslist.ctd[[1]][21:25] <- 5

#make speciation rates and the transition rate optimised
idparsopt.cr <- c(1,3) 
idparsopt.etd <- c(1:5,7) 
idparsopt.ctd <- c(1:5,7)

#fixed dual transition and extinction rate
idparsfix.cr <- c(0,2) 
idparsfix.etd <- c(0,6) 
idparsfix.ctd <- c(0,6) 

#cond set to "proper_cond"
#sampling fraction  = fractions in each state but if don't know how many in each state total use overall sampling fraction for each
sampling_fraction <- rep(length(clusters.j)/305523, c.states) #number of taxa/number of accepted seed plant species 

#get and set starting values
startingpoint <- bd_ML(brts = ape::branching.times(phylo.clusters)) #get starting values with simple birth-death model
intGuessLamba <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0

initparsopt.cr <- c(rep(intGuessLamba, 1), rep(0.25,1))
initparsopt.etd <- c(rep(intGuessLamba, c.states), rep(0.25,1))
initparsopt.ctd <- c(rep(intGuessLamba, c.states), rep(0.25,1))

#set values of dual transition and extinction
parsfix <- c(0,intGuessMu) 

#run models
secsse.cr <- secsse_ml(phylo.clusters, clusters.j, num_concealed_states=5, idparslist.cr, idparsopt.cr, 
                       initparsopt.cr, idparsfix.cr, parsfix, cond="maddison_cond",
                       root_state_weight = "maddison_weights", tol = c(1e-04, 1e-05, 1e-07),
                       sampling_fraction=sampling_fraction, maxiter = 1000 * round((1.25)^length(idparsopt.cr)),
                       methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=T)
saveRDS(secsse.cr, file=paste("Secsse_CR_",duration, "yrs_", unit.type,".Rdata", sep=""))

secsse.etd <- secsse_ml(phylo.clusters, clusters.j, num_concealed_states=5, idparslist.etd, idparsopt.etd, 
                        initparsopt.etd, idparsfix.etd, parsfix, cond="maddison_cond",
                        root_state_weight = "maddison_weights", tol = c(1e-04, 1e-05, 1e-07),
                        sampling_fraction=sampling_fraction, maxiter = 1000 * round((1.25)^length(idparsopt.etd)),
                        methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=T)
saveRDS(secsse.etd, file=paste("Secsse_ETD_",duration, "yrs_", unit.type,".Rdata", sep=""))

secsse.ctd <- secsse_ml(phylo.clusters, clusters.j, num_concealed_states=5, idparslist.ctd, idparsopt.ctd, 
                        initparsopt.ctd, idparsfix.ctd, parsfix, cond="maddison_cond",
                        root_state_weight = "maddison_weights", tol = c(1e-04, 1e-05, 1e-07),
                        sampling_fraction=sampling_fraction, maxiter = 1000 * round((1.25)^length(idparsopt.ctd)),
                        methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=T)
saveRDS(secsse.ctd, file=paste("Secsse_CTD_",duration, "yrs_", unit.type,".Rdata", sep=""))

#calculate AICw
aics <- c(Aic(secsse.cr$ML, k.values[1]),Aic(secsse.etd$ML, k.values[2]),Aic(secsse.ctd$ML, k.values[3]))
delta.aic <- aics-min(aics)
aic.w <- rep(NA,3)
for(l in 1:length(delta.aic)){
  aic.w[l] <- exp(-(delta.aic[l]/2))/sum(exp(-(delta.aic[1]/2)),exp(-(delta.aic[2]/2)),exp(-(delta.aic[3]/2)))
}


#save results
secsse.res.df <- data.frame("model"=character(), "k"=numeric(), "ML"=numeric(), 
                            "AIC"=numeric(), "delta.AIC"=numeric(), "AICw"=numeric(), "dataset"=character(), 
                            "unit.type"=character(), "lambda1"=numeric(), "lambda2"=numeric(),
                            "lambda3"=numeric(),"lambda4"=numeric(),"lambda5"=numeric())

secsse.res.df[1,] <- c("CR", k.values[1], secsse.cr$ML, aics[1], delta.aic[1], aic.w[1], duration, unit.type, secsse.cr$MLpars[[1]][1], rep(NA,4))
secsse.res.df[2,] <- c("ETD", k.values[2], secsse.etd$ML, aics[2], delta.aic[2], aic.w[2], duration, unit.type, secsse.etd$MLpars[[1]][1:5])
secsse.res.df[3,] <- c("CTD", k.values[3], secsse.ctd$ML, aics[3], delta.aic[3],aic.w[3], duration, unit.type, secsse.ctd$MLpars[[1]][c(1,6,11,16,21)])

saveRDS(secsse.res.df, file=paste("SecSSE_cluster_test_results_", duration, "yrs_", unit.type, ".Rdata", sep=""))

### set up comparative data objects for pgls ###
clust.df <- tip.evolution.df[!tip.evolution.df$cluster==0,]
clust.df$cluster <- droplevels(clust.df$cluster)

#change tip names to standardised names to match phylo
clust.df$tip <- mast.df$Species[match(unlist(mast.df$tips_original), clust.df$tip)]
clust.df$tip <- gsub(" ", "_", clust.df$tip)

# for evolutionary metrics
clust.cd <- comparative.data(phylo, clust.df[,c(1,2,5,9)], names.col = tip)

#for masting metrics
mast.df$Species <- gsub(" ", "_", mast.df$Species)
mast.lat.cd <- comparative.data(phylo, mast.df[,c(1,11,12,13,26)], names.col = Species)
mast.spat.cd <- comparative.data(phylo, mast.df[,c(1,16, 26)], names.col = Species)
mast.spat100.cd <- comparative.data(phylo, mast.df[,c(1,15,26)], names.col = Species)

## subset BAMM object with lambda to species in masting for STRAPP analysis ##
mast.clusters.bamm <- mast.df$cluster
names(mast.clusters.bamm) <- mast.df$tips_original
mast.clusters.bamm <- mast.clusters.bamm[names(mast.clusters.bamm) %in% mast.BAMM$tip.label]

## subset BAMM object with beta to species in masting data for STRAPP analysis##
# get the masting data to match the bamm object
seed.clusters <- mast.df$cluster #copy the masting data
names(seed.clusters) <- mast.df$tips_original
seed.clust.BAMM <- subtreeBAMM(seed.BAMM.full, tips = names(seed.clusters), node = NULL)
seed.clusters <- seed.clusters[names(seed.clusters) %in% seed.clust.BAMM$tip.label]

####### model fitting ############
## synchrony - pgls ##
try(s.pgls2 <-  pgls(exp(mean_spat_global)~cluster, data=mast.spat.cd, lambda = "ML"))
try(s100.pgls2 <-  pgls(exp(mean_spat_100)~cluster, data=mast.spat100.cd, lambda = "ML"))

## latitude - pgls ##
try(lat.pgls1 <-  pgls(mean_latitude_absolute~cluster, data=mast.lat.cd, lambda = "ML"))

## lambda (speciation rate) - STRAPP ##
try(strapp.mast.clust <- traitDependentBAMM(mast.BAMM, mast.clusters.bamm, 1000, rate = "speciation", method = "kruskal", logrates = F ,two.tailed = TRUE, nthreads = 1))

## seed size evoul (beta) - STRAPP ##
try(strapp.seed.clust <- traitDependentBAMM(seed.clust.BAMM, seed.clusters, 1000, rate = "trait", method = "kruskal", logrates = F ,two.tailed = TRUE, nthreads = 1))

## DR - pgls ##
try(dr.pgls1 <-  pgls(log(dr)~cluster, data=clust.cd, lambda = "ML"))

### save outputs ###
# macroevoulutionary rates and masting metrics
  #fills in the run info, then model results separately, in case the model didn't work
lamb.row.j <- c("lambda", "STRAPP", rep(NA, 10),length(mast.clusters.bamm))
try(lamb.row.j[7:12] <- c(strapp.mast.clust$p.value, unlist(strapp.mast.clust$estimate)))

seed.row.j <- c("beta", "STRAPP", rep(NA, 10),length(seed.clusters))
try(seed.row.j[7:12] <- c(strapp.seed.clust$p.value, unlist(strapp.seed.clust$estimate)))

dr.row.j <- c("dr", "pgls", rep(NA, 10), nrow(clust.cd$data))
try(dr.row.j[3:7] <- c(summary(dr.pgls1)$fstatistic[2], summary(dr.pgls1)$fstatistic[3],summary(dr.pgls1)$fstatistic[1],
              summary(dr.pgls1)$r.squared, p.value(dr.pgls1)))

s.row.j <- c("s_global", "pgls", rep(NA, 10), nrow(mast.spat.cd$data))
try(s.row.j[3:7] <- c(summary(s.pgls2)$fstatistic[2], summary(s.pgls2)$fstatistic[3],summary(s.pgls2)$fstatistic[1],
             summary(s.pgls2)$r.squared, p.value(s.pgls2)))

s100.row.j <- c("s_100", "pgls", rep(NA, 10), nrow(mast.spat100.cd$data))
try(s100.row.j[3:7] <- c(summary(s100.pgls2)$fstatistic[2], summary(s100.pgls2)$fstatistic[3],summary(s100.pgls2)$fstatistic[1],
                summary(s100.pgls2)$r.squared, p.value(s100.pgls2)))

lat.row.j <- c("latitude", "pgls", rep(NA, 10), nrow(mast.lat.cd$data))
try(lat.row.j[3:7] <- c(summary(lat.pgls1)$fstatistic[2], summary(lat.pgls1)$fstatistic[3],summary(lat.pgls1)$fstatistic[1],
               summary(lat.pgls1)$r.squared, p.value(lat.pgls1)))

 cluster.results <- data.frame("rate"=character(), "method"=character(), "df1"=numeric(), "df2"=numeric(),
                                 "F"=numeric(), "R2"=numeric(), "p"=numeric(), "median_A"=numeric(), "median_B"=numeric(), 
                                "median_C"=numeric(), "median_D"=numeric(), "median_E"=numeric(), "N"=numeric() )
  cluster.results [1,]<- lamb.row.j
  cluster.results [2,]<- seed.row.j
  cluster.results [3,]<- dr.row.j
  cluster.results [4,]<- s.row.j
  cluster.results [5,]<- s100.row.j 
  cluster.results [6,]<- lat.row.j

  cluster.results$dataset <- rep(duration, nrow(cluster.results))
  cluster.results$unit.type <- rep(unit.type, nrow(cluster.results))
  
saveRDS(cluster.results, file=paste("Macroevolutionary_rates_cluster_test_results_", duration, "yrs_", unit.type, ".Rdata", sep=""))

try(rm(s.pgls2))
try(rm(s100.pgls2))
try(rm(lat.pgls1))
try(rm(dr.pgls1))
try(rm(strapp.seed.clust))
try(rm(strapp.mast.clust))

#############################################
#### Fig 2 cluster plots masting metrics ####
#############################################
cluster.cols <- c("#009E73", "#F0E442", "#57B4E9", "#D55E00", "#E69F00")[1:no.clusters]

 plot.name <- paste("Fig 2 cluster plots ", no.clusters, " clusters ", duration, "yrs ", unit.type, ".pdf", sep="")
 pdf(file=plot.name, height=4.5, width=7)

 layout(rbind(c(rep(1,4),rep(2,3)),
              c(rep(1,4),rep(2,3)),
              c(rep(1,4),rep(3,3)),
              c(rep(1,4),rep(3,3))))

 #scatterplot of clusters
 par(mar=c(5.1,4.1,2.1,0.5))
 plot(mast.df$mean_CV~mast.df$mean_AR1, ylab="Mean CVp", xlab="Mean AR1", type="n")
 points(mast.df$mean_CV[mast.df$cluster==1]~mast.df$mean_AR1[mast.df$cluster==1], bg=cluster.cols[1], pch=21)
 points(mast.df$mean_CV[mast.df$cluster==2]~mast.df$mean_AR1[mast.df$cluster==2], bg=cluster.cols[2], pch=21)
 points(mast.df$mean_CV[mast.df$cluster==3]~mast.df$mean_AR1[mast.df$cluster==3], bg=cluster.cols[3], pch=21)
 points(mast.df$mean_CV[mast.df$cluster==4]~mast.df$mean_AR1[mast.df$cluster==4], bg=cluster.cols[4], pch=21)
 points(mast.df$mean_CV[mast.df$cluster==5]~mast.df$mean_AR1[mast.df$cluster==5], bg=cluster.cols[5], pch=21)

 #legend
 text(-0.9, 2.43, "Clusters", pos=4)
 text(rep(-0.85,no.clusters), seq(from=2.35, by=-0.08, length.out = no.clusters), LETTERS[1:no.clusters], pos=4, cex=0.75)
 points(rep(-0.85,no.clusters), seq(from=2.35, by=-0.08, length.out = no.clusters), pch=21, bg=cluster.cols)
 text(-0.92, 2.65, "a)", pos=4)

 # boxplot synchrony
 par(mar=c(4.1,4.1,2.1,1))
 boxplot(mean_spat_global~cluster, data=mast.df, col=cluster.cols, xaxt="n", ylab="Mean synchrony", xlab="")
 axis(1, at=1:no.clusters,labels = LETTERS[1:no.clusters])
 text(0.3, 0.95, "b)", pos=4)
 text(0.3, -0.65, f.fig(s.pgls2), pos=4, cex=0.75)
 text(0.3,-0.80, bquote(paste("R"^"2", "=", .(r.sqr(s.pgls2)))), pos=4, cex=0.75)
 text(0.3,-0.95, paste("p-value =", p.value(s.pgls2)), pos=4, cex=0.75)
 #text(0.3,-0.95, "p-value = 0.90", pos=4, cex=0.75)

 #boxplot mean latitude
 par(mar=c(5.1,4.1,0,1))
 boxplot(mean_latitude_absolute~cluster, data=mast.df, col=cluster.cols, xaxt="n", ylab="Absolute latitude (°)", xlab="Cluster", ylim=c(0,90))
 axis(1, at=1:no.clusters,labels = LETTERS[1:no.clusters])
 text(0.3, 0.95*90, "c)", pos=4)
 text(4, 88, f.fig(lat.pgls1), pos=4, cex=0.75)
 text(4, 81, bquote(paste("R"^"2", "=", .(r.sqr(lat.pgls1)))), pos=4, cex=0.75)
 text(4,74, paste("p-value = ", p.value(lat.pgls1)), pos=4, cex=0.75)

 dev.off()

 #############################################
 #### Fig 5 cluster plots masting metrics ####
 #############################################
fig.name <- paste("Fig 5 macroevolution boxplots by cluster ", no.clusters," clusters ", duration, "yrs ", unit.type, ".pdf", sep="")
 pdf(file=fig.name, height=6, width = 6)
 par(mfrow=c(2,2), mar=c(4,4.1,1,1.5))

 #lambda (speciation rate)
 l <- boxplot(lambda~cluster, data=clust.df, col=cluster.cols, xaxt="n", ylab="Rate of speciation (spp./Ma)", xlab="")
 title(xlab="Cluster", line=2.5)
 axis(1, at=1:no.clusters,labels = LETTERS[1:no.clusters])
 text(0.3, 0.95*max(clust.df$lambda, na.rm = T), "a)", pos=4)
 text(no.clusters-1.75, 2.1, paste("p-value = ", round(strapp.mast.clust$p.value, digits=2)), pos=4, cex=0.75)
 for(i in 1:no.clusters){
   est.number <- grep(i, names(strapp.mast.clust$estimate))
   est <- round(as.numeric(strapp.mast.clust$estimate[est.number]), digits=2)
   if(nchar(est)==3){est <- paste(est, "0", sep="")}
   y.pos <- mean(c(l$stats[3,i],l$stats[4,i]))
   text(i, y.pos, est, cex=0.75)
 }


 # beta (seed mass evolution rate)
 b <- boxplot(beta~cluster, data=clust.df, col=cluster.cols, log="y", xaxt="n", ylab="Rate of seed mass change", xlab="")
 axis(1, at=1:no.clusters,labels = LETTERS[1:no.clusters])
 title(xlab="Cluster", line=2.5)
 text(0.3, 0.95*max(clust.df$beta, na.rm = T), "b)", pos=4)
 text(no.clusters-1.7, 0.009, paste("p-value = ", round(strapp.seed.clust$p.value, digits=2)), pos=4, cex=0.75)
 for(i in 1:no.clusters){
   est.number <- grep(i, names(strapp.seed.clust$estimate))
   est <- round(as.numeric(strapp.seed.clust$estimate[est.number]), digits=2)
   if(nchar(est)==3){est <- paste(est, "0", sep="")}
   y.pos <- mean(c(b$stats[3,i],b$stats[4,i]))
   text(i, y.pos, est, cex=0.75)
 }
 # DR (diversification rate)
 boxplot(dr~cluster, data=clust.df, col=cluster.cols, log="y",  xaxt="n", ylab="Diversification rate metric (DR, spp./Ma)",
         xlab="", ylim=c(00.008,max(clust.df$dr)))
 axis(1, at=1:no.clusters,labels = LETTERS[1:no.clusters])
 title(xlab="Cluster", line=2.5)
 text(0.3, 0.95*max(clust.df$dr), "c)", pos=4)
 text(no.clusters-4.7, exp(-4.2), f.fig(dr.pgls1), pos=4, cex=0.75)
 text(no.clusters-4.7,exp(-4.55), bquote(paste("R"^"2", "=", .(r.sqr(dr.pgls1)))), pos=4, cex=0.75)
 text(no.clusters-4.7,exp(-4.9), paste("p-value = ", p.value(dr.pgls1)), pos=4, cex=0.75)

 dev.off()

}#end of j loop

#### combine dataframes ####
k.result.dfs <-list.files(pattern="Phylogenetic_signal_results_.")
cor.result.dfs <- list.files(pattern="Phylogenetic_correlation_results_.")
trait.divn.dfs <- list.files(pattern="Trait_dependent_diversification_results_masting_.")[-c(1,7,13,19, 25:27)]
quasse.output.dfs <-list.files(pattern="Quasse_results_masting_.") 
bias.dfs <- list.files(pattern="Macroevolutionary_bias_mastree_data_.")
cluster.res.dfs <- list.files(pattern="Macroevolutionary_rates_cluster_test_results")
secsse.res.dfs <- list.files(pattern="SecSSE_cluster_test_results_.")


for(i in 1:length(secsse.res.dfs)){
  k.res.i <- readRDS(k.result.dfs[i])
  cor.res.i <- readRDS(cor.result.dfs[i])
  trait.divn.i <- readRDS(trait.divn.dfs[i])
  quasse.output.i <- readRDS(quasse.output.dfs[i])
  bias.i <- readRDS(bias.dfs[i])
  cluster.res.i <- readRDS(cluster.res.dfs[i])
  secsse.res.i <- readRDS(secsse.res.dfs[i])
  if(i ==1){
    k.res.master <- k.res.i
    cor.res.master <- cor.res.i
    trait.divn.master <- trait.divn.i
    quasse.output.master <- quasse.output.i
    bias.master <- bias.i
    cluster.res.master <- cluster.res.i
    secsse.res.master <- secsse.res.i
  }else{
    k.res.master <- rbind(k.res.master, k.res.i)
    cor.res.master <- rbind(cor.res.master, cor.res.i)
    trait.divn.master <- rbind(trait.divn.master, trait.divn.i)
    quasse.output.master <- rbind(quasse.output.master, quasse.output.i)
    bias.master <- rbind(bias.master,bias.i)
    cluster.res.master <- rbind(cluster.res.master, cluster.res.i)
    secsse.res.master <- rbind(secsse.res.master, secsse.res.i)
  }
}
saveRDS(k.res.master, file="Phylogenetic_signal_resultsMaster.Rdata")
saveRDS(cor.res.master, file="Phylogenetic_correlation_resultsMaster.Rdata")
saveRDS(trait.divn.master, file="Trait_dependent_diversification_results_mastingMaster.Rdata")
saveRDS(quasse.output.master, file="Quasse_results_mastingMaster.Rdata")
saveRDS(bias.master, file="Macroevolutionary_bias_mastree_dataMaster.Rdata")
saveRDS(cluster.res.master, file="Macroevolutionary_rates_cluster_test_resultsMaster.Rdata")
saveRDS(secsse.res.master, file="SecSSE_cluster_test_resultsMaster.Rdata")

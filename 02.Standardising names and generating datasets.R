## checking names in MASTREE+ data, Smith and Brown GBOTB phylogeny, Smith and Brown BAMM object, and seed BAMM object ##
# does name check and filters data to units of interest and one unit per species for 3:6 year duration thresholds and makes corresponding phylgeny
# Written by Esther Dale, Feb 2021

library(Taxonstand)
library(taxize)
library(dplyr)
library(ape)

# read in masting data
load("Summary_All_indices_3_years_pooled.RData")
mast.df <-sum_all_indices_pooled
mast.names.raw <- unique(mast.df$Species)
mast.names.raw <- gsub("Stenosmum acreanum", "Stenostomum acreanum", mast.names.raw) #this species seems to have trouble matching and drop out of gnrs list
mast.names.raw <- gsub("Garcinia eugeniaefolia", "Garcinia eugeniifolia", mast.names.raw) #this species doesn't match on TPL or GNRS

# run through TPL
mast.names.checked <- TPL(splist = mast.names.raw)

#run taxa that didn't match on TPL (e.g. mispelled genera) through gnrs
mast.names.rechecked <- gnr_resolve(sci=as.character(mast.names.checked$Taxon[mast.names.checked$Plant.Name.Index==FALSE]), best_match_only = T)

#then back through TPL
mast.names.rechecked2 <- TPL(splist = mast.names.rechecked$matched_name)
mast.names.rechecked2 <- left_join(mast.names.rechecked2, mast.names.rechecked[,c(1,3)], by=c("Taxon"="matched_name"))
mast.names.rechecked2$Taxon <-mast.names.rechecked2$user_supplied_name #change back to original name so spelling change is recorded
mast.names.checked<- rbind(mast.names.checked[mast.names.checked$Plant.Name.Index==TRUE,], mast.names.rechecked2[,1:26])

#filter out unsuitable taxa
mast.names.accepted <- mast.names.checked[mast.names.checked$Tax_res=="Species",]
mast.names.accepted <- mast.names.accepted[complete.cases(mast.names.accepted$New.Species),]
saveRDS(mast.names.accepted, file="Mast names accepted full.Rdata")

#read in phylogeny
tree.tr <- read.tree("GBOTB.tre")
  #phylogeny from Smith & Brown 2018

tip.labs.df <- as.data.frame(tree.tr$tip.label, stringsAsFactors = F)
colnames(tip.labs.df) <- "tip_lab"
tip.labs.df$genus <- sapply(strsplit(tip.labs.df$tip_lab, split = "_"), `[`, 1) #make genus column
tip.labs.full <- tip.labs.df

#subset to genera in masting data
tip.labs.df <- tip.labs.df[tip.labs.df$genus %in% unique(c(mast.names.accepted$Genus,mast.names.accepted$New.Genus)),]
tip.labs.df$tip_lab <- gsub("_", " ", tip.labs.df$tip_lab) #remove underscore

# run through TPL
  # done in sections so progress isn't lost with connectivity issues
index.list <- list(1:1000, 1001:2000, 2001:3000, 3001:4000, 4001:5000, 5001:6000, 6001:7000, 7001:nrow(tip.labs.df))
 for(i in 1:length(index.list)){
  tip.labs.partial <- TPL(splist = tip.labs.df$tip_lab[index.list[[i]]])
if(i==1){tip.labs.checked <- tip.labs.partial}else{
  tip.labs.checked <- rbind(tip.labs.checked,tip.labs.partial)
} #if else
} #i loop

saveRDS(tip.labs.checked, file="TPL output GBOTB names masting genera.Rdata")

#run taxa that didn't match on TPL (e.g. mispelled genera) through gnr
tip.labs.rechecked <- gnr_resolve(sci=as.character(tip.labs.checked$Taxon[tip.labs.checked$Plant.Name.Index==FALSE]), best_match_only = T)
  #remove species that are genus only
tip.labs.rechecked$genus <- sapply(tip.labs.rechecked$Name_matched, function (x) length(unlist(strsplit(as.character(x), "\\W+"))))
tip.labs.rechecked <- tip.labs.rechecked[tip.labs.rechecked$genus>1,]

#then back through TPL
tip.labs.rechecked2 <- TPL(splist = tip.labs.rechecked$Name_matched)
tip.labs.rechecked2 <- left_join(tip.labs.rechecked2, tip.labs.rechecked[,c(1,2)], by=c("Taxon"="Name_matched"))
tip.labs.rechecked2$Taxon <- tip.labs.rechecked2$Name_submitted #change back to original name so spelling change is recorded
tip.labs.checked<- rbind(tip.labs.checked[tip.labs.checked$Plant.Name.Index==TRUE,], tip.labs.rechecked2[,1:26])
  #check the species that have multiple possible corrections if they're in the masting data

#match masting data names to phylo names
mast.names.accepted$binomial <- paste(mast.names.accepted$New.Genus, mast.names.accepted$New.Species, sep=" ")
tip.labs.checked$binomial <- paste(tip.labs.checked$New.Genus, tip.labs.checked$New.Species, sep=" ")

no.match <- mast.names.accepted[which(!mast.names.accepted$binomial %in% tip.labs.checked$binomial),]
no.match$genus.change <- (!no.match$Genus==no.match$New.Genus)
no.match$binomial[no.match$genus.change==T] %in% tip.labs.full$tip_lab #see if any species that had a genus change could have been missed

#subset to species that are in masting data and phylo
mast.names.accepted <- mast.names.accepted[mast.names.accepted$binomial %in% tip.labs.checked$binomial,]
mast.names.accepted <- mast.names.accepted[(!mast.names.accepted$Taxon=="Acacia craspedocarpa x aneura"),] #remove a hybrid
  
#subset masting data to matched species and change to accepted names
mast.df <- mast.df[mast.df$Species %in% mast.names.accepted$Taxon,]
mast.df$Species.original <- mast.df$Species

for(i in 1:nrow(mast.names.accepted)){
  name.i <- mast.names.accepted$Taxon[i] #ith original name
  mast.df$Species <- gsub(name.i, mast.names.accepted$binomial[i], mast.df$Species)
}

### prune phylogeny to masting species
tip.labs.checked <- tip.labs.checked[tip.labs.checked$binomial %in% mast.names.accepted$binomial,] #subset to species in masting data

# filter out subspecies that aren't needed
subspecies.list <- unique(tip.labs.checked$binomial[!tip.labs.checked$Infraspecific==""]) #list of binomials that have subspecies
subspecies.to.keep <- subspecies.list[!subspecies.list %in% unique(tip.labs.checked$binomial[tip.labs.checked$Infraspecific==""])] #remove subspecies that have a species-level tip
subspecies.to.keep <- tip.labs.checked$Taxon[tip.labs.checked$binomial %in% subspecies.to.keep]
subspecies.to.keep <- subspecies.to.keep[7] #select which ones to keep
tip.labs.checked <- tip.labs.checked[tip.labs.checked$Infraspecific==""|tip.labs.checked$Taxon%in% subspecies.to.keep,]

#filter out duplicates
duplicates.list <- unique(tip.labs.checked$binomial[which(duplicated(tip.labs.checked$binomial))])
duplicates.to.keep <- tip.labs.checked$Taxon[tip.labs.checked$binomial%in% duplicates.list & tip.labs.checked$Species== tip.labs.checked$New.Species]
duplicates.rows <- which(tip.labs.checked$binomial %in% duplicates.list & (!tip.labs.checked$Taxon %in% duplicates.to.keep))
tip.labs.checked <- tip.labs.checked[(!1:nrow(tip.labs.checked) %in% duplicates.rows),]

## filter data to variables and units of interest ##
mast.names.accepted <- mast.names.accepted[mast.names.accepted$binomial %in% tip.labs.checked$binomial,]
mast.df <- mast.df[mast.df$Species %in% mast.names.accepted$binomial,]
#filter out variables we're not interested in
mast.df <- mast.df[!mast.df$Unit %in% c("% aboveground biomass", "% branch fruiting", "% fruit score",  "% indivs reproducing", "1000kcal/ha",  "cones/branch", "cones/shoot", "fruits/branch", "g/branch", "index", "kg", "kJ/tree", "m3","MJ/ha", "seeds", "seeds/branch"),]
#filter data to one unit type per species
  #picks the variable with the most time series (n_timeseries)
species <- unique(mast.df$Species)
mast.df$keep <- rep(0, nrow(mast.df))
for(i in 1:length(species)){
  rows.i <- which(mast.df$Species==species[i])
  if(length(rows.i>1)){
    to.keep.i <- rows.i[which(mast.df$n_timeseries[rows.i]==sort(mast.df$n_timeseries[rows.i], decreasing = T)[1])[1]]  #which row has highest count
    mast.df$keep[to.keep.i] <- 1
  }else{
    mast.df$keep[rows.i] <- 1
  }
}
mast.df <- mast.df[mast.df$keep==1,]

write.csv(mast.names.accepted,file="Masting species accepted names.csv")

saveRDS(mast.df[,1:22], file="Full masting data names checked.Rdata")

# check names of bamm object in sections
library(BAMMtools)
library(matrixStats)
Smith.Brown.BAMM <- readRDS("SmithBrown_BAMM_vascular.RDS")
#BAMM object from Igea & Tanentzap 2020 https://github.com/javierigea/LDGplants_rates
bamm.labs.df <- as.data.frame(Smith.Brown.BAMM$tip.label, stringsAsFactors = F)
colnames(bamm.labs.df) <- "tip_lab"
bamm.labs.df$genus <- sapply(strsplit(bamm.labs.df$tip_lab, split = "_"), `[`, 1) #make genus column
bamm.labs.df$tip_lab <- gsub("_", " ", bamm.labs.df$tip_lab)
bamm.labs.full <- bamm.labs.df

index.list <- list(1:1000, 1001:2000, 2001:3000, 3001:4000, 4001:5000, 5001:6000, 6001:7000, 7001:73934)
for(i in 1:length(index.list)){
  bamm.labs.partial <- TPL(splist = bamm.labs.df$tip_lab[index.list[[i]]])
  if(i==1){bamm.labs.checked <- bamm.labs.partial}else{
    bamm.labs.checked <- rbind(bamm.labs.checked,bamm.labs.partial)
  } #if else
} #i loop
saveRDS(bamm.labs.checked, file="TPL output Smith and Brown BAMM names masting genera.Rdata")

bamm.labs.checked$binomial <- paste(bamm.labs.checked$New.Genus, bamm.labs.checked$New.Species, sep=" ")

#prune tree to tip labs checked
tip.labs.checked$tip_lab <- gsub(" ", "_", tip.labs.checked$Taxon) #add back in underscore
tree.tr <- keep.tip(tree.tr, tip=tip.labs.checked$tip_lab)
#change tip labels to TPL checked binomials
tree.tr$tip.label <- tip.labs.checked$binomial[match(tip.labs.checked$tip_lab, tree.tr$tip.label)]
#save phylo
write.tree(tree.tr, file="Masting phylogeny.nwk")

# read in seed data
seed.data <- readRDS("Trait_edata_50shifts_1000samples.Rsave")
#BAMM object from Igea et al. 2017 https://github.com/javierigea/seed_size
seed.names <-as.data.frame(seed.data$tip.label)
colnames(seed.names) <- "tip.label"
seed.names$tip.label <- as.character(seed.names$tip.label)
seed.names$tip.label <- gsub("_", " ", seed.names$tip.label)

#check names with TPL
seed.names.checked<- TPL(splist =seed.names$tip.label[seed.names$genus %in% mast.names.accepted$Genus])

#subset to species in accepted masting species list
seed.names.checked$binomial <- paste(seed.names.checked$New.Genus, seed.names.checked$New.Species, sep=" ")
seed.names.checked <- seed.names.checked[seed.names.checked$binomial %in% mast.names.accepted$binomial,]
saveRDS(seed.names.checked, file="Seed names checked.Rdata")

#### check names and make name checked phylogenies for the different duration datasets####
# load data #
mast3.df <- mast.df[,1:22] #remove keep column to match others
load("Summary_All_indices_4_years_pooled.RData")
mast4.df <- sum_all_indices_pooled
load("Summary_All_indices_5_years_pooled.RData")
mast5.df <- sum_all_indices_pooled
load("Summary_All_indices_6_years_pooled.RData")
mast6.df <- sum_all_indices_pooled
df.list <- list(mast3.df,mast4.df,mast5.df, mast6.df)

for(i in 1:length(df.list)){
  df.i <- df.list[[i]]
  df.i <- as.data.frame(df.i)
  df.i$Species.original <- df.i$Species
  #find checked name
  df.i$Species  <- mast.df$Species[match(unlist(df.i$Species.original), mast.df$Species.original)]
  #subset to species that have a checked name
  df.i <- df.i[complete.cases(df.i$Species),]
  
  #filter out variables we're not interested in
  df.i <- df.i[!df.i$Unit %in% c("% aboveground biomass", "% branch fruiting", "% fruit score",  "% indivs reproducing", "1000kcal/ha",  "cones/branch", "cones/shoot", "fruits/branch", "g/branch", "index", "kg", "kJ/tree", "m3","MJ/ha", "seeds", "seeds/branch"),]
  df.i$keep <- rep(0, nrow(df.i))
  
  # make version of df with are and individual based units
  df.i.area <- df.i[(df.i$Unit %in% c("cones/m2", "cones/plot", "seed/m2BA", "fruits/m2", "g/crown area m2", "g/m2","g/plot","seeds/crown area m2", "seeds/m2", "seeds/plot")),]
  df.i.indiv <- df.i[(df.i$Unit %in% c("cones/individual", "g/individual", "seeds/individual")),]
  
  #make version with mass and counts based units
  df.i.mass <- df.i[(df.i$Unit %in% c("g/crown area m2", "g/m2","g/plot", "g/individual")),]
  df.i.count <- df.i[(df.i$Unit %in% c("cones/m2", "cones/plot", "seed/m2BA", "fruits/m2", "seeds/crown area m2", "seeds/m2", "seeds/plot", "cones/individual", "seeds/individual")),]

  
   #filter data to one unit type per species
    #picks the variable with the most time series (n_timeseries)
  species <- unique(df.i$Species)
  for(j in 1:length(species)){
    rows.j <- which(df.i$Species==species[j])
    rows.j.area <- which(df.i.area$Species==species[j])
    rows.j.indiv <- which(df.i.indiv$Species==species[j])
    rows.j.mass <- which(df.i.mass$Species==species[j])
    rows.j.count <- which(df.i.count$Species==species[j])
    
    #filter for all units
    if(length(rows.j>1)){
      to.keep.j <- rows.j[which(df.i$n_timeseries[rows.j]==sort(df.i$n_timeseries[rows.j], decreasing = T)[1])[1]]  #which row has highest count
      df.i$keep[to.keep.j] <- 1
    }else{
      df.i$keep[rows.j] <- 1
    }
    #filter for area units
    if(length(rows.j.area>1)){
      to.keep.j.area <- rows.j.area[which(df.i.area$n_timeseries[rows.j.area]==sort(df.i.area$n_timeseries[rows.j.area], decreasing = T)[1])[1]]  #which row has highest count
      df.i.area$keep[to.keep.j.area] <- 1
    }else{
      df.i.area$keep[rows.j.area] <- 1
    }
    #filter for indiv units
    if(length(rows.j.indiv>1)){
      to.keep.j.indiv <- rows.j.indiv[which(df.i.indiv$n_timeseries[rows.j.indiv]==sort(df.i.indiv$n_timeseries[rows.j.indiv], decreasing = T)[1])[1]]  #which row has highest count
      df.i.indiv$keep[to.keep.j.indiv] <- 1
    }else{
      df.i.indiv$keep[rows.j.indiv] <- 1
    }
    #filter for mass units
    if(length(rows.j.mass>1)){
      to.keep.j.mass <- rows.j.mass[which(df.i.mass$n_timeseries[rows.j.mass]==sort(df.i.mass$n_timeseries[rows.j.mass], decreasing = T)[1])[1]]  #which row has highest count
      df.i.mass$keep[to.keep.j.mass] <- 1
    }else{
      df.i.mass$keep[rows.j.mass] <- 1
    }
    #filter for counts units
    if(length(rows.j.count>1)){
      to.keep.j.count <- rows.j.count[which(df.i.count$n_timeseries[rows.j.count]==sort(df.i.count$n_timeseries[rows.j.count], decreasing = T)[1])[1]]  #which row has highest count
      df.i.count$keep[to.keep.j.count] <- 1
    }else{
      df.i.count$keep[rows.j.count] <- 1
    }
  }
  df.i <- df.i[df.i$keep==1,]
  df.i.area <- df.i.area[df.i.area$keep==1,]
  df.i.indiv <- df.i.indiv[df.i.indiv$keep==1,]
  df.i.mass <- df.i.mass[df.i.mass$keep==1,]
  df.i.count <- df.i.count[df.i.count$keep==1,]
  
  #save
  df.name <- paste("Masting data names checked", (3:6)[[i]], "yrs.Rdata", sep=" ")
  df.name.area <- paste("Masting data names checked", (3:6)[[i]], "yrs area units only.Rdata", sep=" ")
  df.name.indiv <- paste("Masting data names checked", (3:6)[[i]], "yrs individual units only.Rdata", sep=" ")
  df.name.mass <- paste("Masting data names checked", (3:6)[[i]], "yrs mass units only.Rdata", sep=" ")
  df.name.count <- paste("Masting data names checked", (3:6)[[i]], "yrs count units only.Rdata", sep=" ")
  
   saveRDS(df.i, df.name)
  saveRDS(df.i.area, df.name.area)
  saveRDS(df.i.indiv, df.name.indiv)
  saveRDS(df.i.mass, df.name.mass)
  saveRDS(df.i.count, df.name.count)
  
  #subset phylogeny to species in df.i
  tree.i <- keep.tip(tree.tr, tip = df.i$Species)
  tree.name <- paste("Masting phylogeny ", (3:6)[i], "yrs.nwk", sep="")
  write.tree(tree.i, file=tree.name)
}

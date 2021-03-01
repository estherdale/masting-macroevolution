# species ages function from Tanentzap et al.(2020). Does evolutionary history correlate with contemporary extinction risk by influencing range size dynamics?. The American Naturalist, 195(3), 569-576.

calculate_tip_ages<-function(tree){
  node.age(tree)->phy.age
  cbind(phy.age$edge,phy.age$age, tree$edge.length)->BL.position
  max(phy.age$age)-BL.position[,3]->dist.tip
  cbind(BL.position,dist.tip)->BL.positions
  BL.positions[,5]+BL.positions[,4]->ages
  cbind(BL.positions,ages)->BL.positions
  as.data.frame(BL.positions)->node.ages
  names(node.ages)<-c("parental.node","daughter.node","dist.root","BL","dist.tip","mrca.age")
  ## node.ages is a data frame listing as variables the identity of parental and
  #daughter nodes, the distance from the root and from the present of each node,
  #the branch length and the age of the most recent common ancestor
  node.ages[node.ages[,2]<length(tree$tip)+1,]->species.ages
  row.names(species.ages)<-tree$tip[species.ages$daughter.node]
  ## species ages is node.ages data frame reduced to the tips (species)
  species.ages<-species.ages[order(row.names(species.ages)),]
  output.table<-as.data.frame(cbind(row.names(species.ages),species.ages$mrca.age))
  colnames(output.table)<-c('tip','tip.age')
  return(output.table)
}
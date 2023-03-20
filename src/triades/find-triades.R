####################FIND TRIADES####################
siblings <- edges[edges$sameParent & !edges$sameStrain & edges$seq_shared_best_traj >= 0,]
families <- siblings%>%
  group_by(commonParent)%>%
  summarise(IBDmax = max(seq_shared_best_traj), IBDmin = min(seq_shared_best_traj))
families <- families[families$IBDmin < 0.2 & families$IBDmax >= parentalmin,]

allchildren <- data.frame("ID" = NA, "parent1" = NA, "parent2" = NA, 
                          "parentsIBD" = NA, "misbases" = NA,
                          "parent1IBD" = NA, "parent2IBD" = NA, 
                          "parent1bases" = NA, "parent2bases" = NA, 
                          "samebases" = NA, "family" = NA)

for (family in families$commonParent){
  siblings.family <- siblings[siblings$commonParent==family & siblings$seq_shared_best_traj < parentalmax,]
  siblings.family.unrelated <- siblings.family[siblings.family$seq_shared_best_traj < different,]
  siblings.family.unrelated[,c("ID1", "ID2")] <- data.frame(pmin(siblings.family.unrelated$ID1, 
                                                              siblings.family.unrelated$ID2),
                                                         pmax(siblings.family.unrelated$ID1,
                                                              siblings.family.unrelated$ID2)) 
  
  siblings.family.related <- siblings.family[siblings.family$seq_shared_best_traj > parentalmin,]
  
  wayone <- siblings.family.related%>%
    group_by(ID1)%>%
    summarise(expand.grid(ID2, ID2, stringsAsFactors = F))
  waytwo <- siblings.family.related%>%
    group_by(ID2)%>%
    summarise(expand.grid(ID1, ID1, stringsAsFactors = F))
  
  colnames(wayone) <- c("ID", "parent1", "parent2")
  colnames(waytwo) <- c("ID", "parent1", "parent2")
  
  children <- rbind(wayone, waytwo)
  
  children[,c("parent1", "parent2")] <- data.frame(pmin(children$parent1, 
                                                        children$parent2),
                                                   pmax(children$parent1,
                                                        children$parent2))
  children$parentsIBD <- apply(children, 1,
                        FUN = function(x){
                          value <- siblings.family.unrelated$seq_shared_best_traj[siblings.family.unrelated$ID1 == x[2] & siblings.family.unrelated$ID2 == x[3]]
                          if (length(value) > 0){
                            return(value)
                          }else{
                            return(NA)
                          }
                        })
  children <- children[!is.na(children$parentsIBD),]
  children <- children[!duplicated(children),]
  children$misbases <- apply(children,1,
                                FUN = function(x){
                                  dfunkown <- parentsdiff(chrombarcodes, x[1], x[2], x[3])
                                  return(sum(dfunkown$end - dfunkown$start))
                                })
  children <- children[order(children$ID, children$misbases),]
  children.best <- children[!duplicated(children$ID),]
  children.best[order(children.best$parent1, children.best$parent2),]
  children.best$parent1IBD <- apply(children.best, 1,
                               FUN = function(x){
                                 value <- siblings.family$seq_shared_best_traj[siblings.family$ID1 == x[1] & siblings.family$ID2 == x[2]]
                                 value2 <- siblings.family$seq_shared_best_traj[siblings.family$ID1 == x[2] & siblings.family$ID2 == x[1]]
                                 if (length(value) == 1){
                                   return(value)
                                 }else{
                                   return(value2)
                                 }
                               })
  children.best$parent2IBD <- apply(children.best, 1,
                               FUN = function(x){
                                 value <- siblings.family$seq_shared_best_traj[siblings.family$ID1 == x[1] & siblings.family$ID2 == x[3]]
                                 value2 <- siblings.family$seq_shared_best_traj[siblings.family$ID1 == x[3] & siblings.family$ID2 == x[1]]
                                 if (length(value) == 1){
                                   return(value)
                                 }else{
                                   return(value2)
                                 }
                               })
  children.best$parent1bases <- apply(children.best,1,
                                      FUN = function(x){
                                        subdata <- subduo(chrombarcodes, x[1], x[2])
                                        return(sum(subdata$end[subdata$different==0] - subdata$start[subdata$different==0]))
                                      })
  children.best$parent2bases <- apply(children.best,1,
                                      FUN = function(x){
                                        subdata <- subduo(chrombarcodes, x[1], x[3])
                                        return(sum(subdata$end[subdata$different==0] - subdata$start[subdata$different==0]))
                                      })
  children.best$samebases <- apply(children.best,1,
                               FUN = function(x){
                                 dfsame <- parentssame(chrombarcodes, x[1], x[2], x[3])
                                 return(sum(dfsame$end - dfsame$start))
                               })
  children.best$family <- family
  allchildren <- rbind(allchildren, children.best)
}

allchildren <- allchildren[!is.na(allchildren$ID),]
########################################

########################SAVE DATA###################
write.table(allchildren, "out/triades.tsv", quote = F, row.names = F)
########################################

########################PLOT###################

for (i in 1:nrow(allchildren)){
  child <- as.character(allchildren[i,1])
  parent1 <- as.character(allchildren[i,2])
  parent2 <- as.character(allchildren[i,3])
  gg <- parentsmapv3(chrombarcodes, 
             child, 
             parent1, parent2,
             chlen = chromlength)
  png(paste0("out/offsprings/", nodes$strain[nodes$ID==child], "/",
             child, "_from_", parent2, "-", parent1, ".png"), 
      width = 8, height = 5, res = 200, units = "in")
  plot(gg)
  dev.off()
  gg <- parentsmapv3(chrombarcodes, 
                     child, 
                     parent2, parent1,
                     chlen = chromlength)
  png(paste0("out/offsprings/", nodes$strain[nodes$ID==child], "/",
             child, "_from_", parent1, "-", parent2, ".png"), 
      width = 8, height = 5, res = 200, units = "in")
  plot(gg)
  dev.off()
  allnextchildren <- allchildren[-c(1:i),1]
  for (nextchild in allnextchildren){
    gg <- parentsmapv3(chrombarcodes, 
                       child, 
                       nextchild, nextchild,
                       chlen = chromlength, mode = 1)
    gg
    ggsave(paste0("out/siblings/", child, "_with_", nextchild, ".png"), 
           plot = gg, width = 8, height = 5)
  }
}

########################################
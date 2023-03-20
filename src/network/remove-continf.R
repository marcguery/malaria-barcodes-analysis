###################DATA###################
#get clusters within individuals from network-doi-strains
continf <- nodesfilteredclusters[nodesfilteredclusters$IBD==0.9,]
continf.summary <- continf %>%
  group_by(individual, cluster) %>%
  summarise(earlierind = min(ID), numind = length(ID))
continf.summary <- continf.summary[continf.summary$cluster != "00" & continf.summary$numind>1,]
continf <- merge(continf, continf.summary)
edges.nocontinf <- merge(edges, continf[,c("ID", "earlierind")],
                         by.x = "ID1", by.y = "ID", all.x = TRUE)
edges.nocontinf <- merge(edges.nocontinf, continf[,c("ID", "earlierind")],
                                 by.x = "ID2", by.y = "ID", all.x = TRUE, 
                         suffixes = c("1","2"))
edges.nocontinf$earlierind1[is.na(edges.nocontinf$earlierind1)] <- edges.nocontinf$ID1[is.na(edges.nocontinf$earlierind1)]
edges.nocontinf$earlierind2[is.na(edges.nocontinf$earlierind2)] <- edges.nocontinf$ID2[is.na(edges.nocontinf$earlierind2)]
edges.nocontinf$Continf <- edges.nocontinf$earlierind1!=edges.nocontinf$ID1 | edges.nocontinf$earlierind2!=edges.nocontinf$ID2
nodes.nocontinf <- nodes[nodes$ID%in%c(edges.nocontinf$ID1, edges.nocontinf$ID2),]

write.table(edges.nocontinf,
            file = "out/edges-nocontinf.tsv",
            quote = F, row.names = F)
write.table(nodes.nocontinf,
            file = "out/nodes-nocontinf.tsv",
            quote = F, row.names = F)

edgesfiltered.nocontinf <- merge(edgesfiltered, continf[,c("ID", "earlierind")],
                                 by.x = "ID1", by.y = "ID", all.x = TRUE)
edgesfiltered.nocontinf <- merge(edgesfiltered.nocontinf, continf[,c("ID", "earlierind")],
                                 by.x = "ID2", by.y = "ID", all.x = TRUE, suffixes = c("1","2"))
edgesfiltered.nocontinf$earlierind1[is.na(edgesfiltered.nocontinf$earlierind1)] <- edgesfiltered.nocontinf$ID1[is.na(edgesfiltered.nocontinf$earlierind1)]
edgesfiltered.nocontinf$earlierind2[is.na(edgesfiltered.nocontinf$earlierind2)] <- edgesfiltered.nocontinf$ID2[is.na(edgesfiltered.nocontinf$earlierind2)]
edgesfiltered.nocontinf$Continf <- edgesfiltered.nocontinf$earlierind1!=edgesfiltered.nocontinf$ID1 | edgesfiltered.nocontinf$earlierind2!=edgesfiltered.nocontinf$ID2
nodesfiltered.nocontinf <- nodesfiltered[nodesfiltered$ID%in%c(edgesfiltered.nocontinf$ID1, edgesfiltered.nocontinf$ID2),]

write.table(edgesfiltered.nocontinf,
            file = "out/filtered-edges-nocontinf.tsv",
            quote = F, row.names = F)
write.table(nodesfiltered.nocontinf,
            file = "out/filtered-nodes-nocontinf.tsv",
            quote = F, row.names = F)

relatednessscore1 <- edges.nocontinf[edges.nocontinf$Continf==FALSE,]%>%
  group_by(ID1, )%>%
  summarize(score1 = max(fract_sites_IBD))
relatednessscore2 <- edges.nocontinf[edges.nocontinf$Continf==FALSE,]%>%
  group_by(ID2)%>%
  summarize(score2 = max(fract_sites_IBD))
relatednessscore <- merge(relatednessscore1, relatednessscore2, 
                          by.x = "ID1",by.y = "ID2", 
                          all = T)
relatednessscore$score <- pmax(relatednessscore$score1, 
                               relatednessscore$score2, na.rm = T)
relatednessscore <- merge(relatednessscore,  nodes, by.x = "ID1", by.y = "ID")


nrow(relatednessscore)
length(which(relatednessscore$score < 0.5))/nrow(relatednessscore)
length(which(relatednessscore$score < 0.9))/nrow(relatednessscore)
summary(relatednessscore$score)

relatednessscore.date <- relatednessscore%>%
  group_by(date)%>%
  summarize(numobs = length(score),
            propshared = length(which(score > 0.9))/length(score))

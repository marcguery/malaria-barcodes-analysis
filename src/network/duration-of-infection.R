#Keep only the edgesfiltered connecting the same individual or 2 individuals belonging to this dataset
edgesfiltered.lasting <- edgesfiltered[edgesfiltered$sameCluster & edgesfiltered$sameIndividual & edgesfiltered$commonCluster!="00",]

nodesfiltered.lasting <- edgesfiltered.lasting %>%
  group_by(commonIndividual, commonCluster)%>%
  summarise(days=max(dayElapsed),
            start = min(pmin(date1, date2)),
            end = max(pmax(date1, date2)))
colnames(nodesfiltered.lasting) <- c("individual", "cluster", "days", "start", "end")

edgesfiltered.subsist <- edgesfiltered[edgesfiltered$sameCluster & edgesfiltered$commonCluster!="00",]
edgesfiltered.subsist <- edgesfiltered.subsist %>%
  group_by(commonCluster)%>%
  summarise(days = max(dayElapsed),
            start = min(pmin(date1, date2)),
            end=max(pmax(date1, date2)))
edgesfiltered.subsist <- edgesfiltered.subsist[edgesfiltered.subsist$days > 0,]
nodesfiltered.subsist <- nodesfiltered[nodesfiltered$cluster%in%edgesfiltered.subsist$commonCluster,]

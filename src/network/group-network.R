#######################COHORT REMOVAL########################
#We remove cohort to remove
# effect of picking individuals
# Though we keep members of 1612 since the cohort was made based on 
# infected individuals of 1612 from a randomly sampled population
edgesfiltered.all.nocoh <- edges[edges$date1 < "2017-01-01" & edges$date2 < "2017-01-01" & edges$fract_sites_IBD > -1,]
edgesfiltered.nocoh <- edgesfiltered[edgesfiltered$date1 < "2017-01-01" & edgesfiltered$date2 < "2017-01-01",]
nodesfiltered.all.nocoh <- nodes[nodes$date < "2017-01-01",]
seasons.nocoh <- seasons[seasons$date1 < "2017-01-01" & seasons$date2 < "2017-01-01",]
seasons.combination.nocoh <- seasons.combination[seasons.combination$date1start < "2017-01-01" &
                                                   seasons.combination$date1end < "2017-01-01" &
                                                   seasons.combination$date2start < "2017-01-01" &
                                                   seasons.combination$date2end < "2017-01-01",]
##############################################

edgesfiltered.all.nocoh <- edgesfiltered.all.nocoh[,-which(colnames(edgesfiltered.all.nocoh)=="fract_sites_IBD")]
edgesfiltered.all.nocoh <- merge(edgesfiltered.all.nocoh, edgesfiltered.nocoh[,c("ID1", "ID2", "fract_sites_IBD")],
                           all.x = TRUE)
edgesfiltered.all.nocoh$fract_sites_IBD[!is.na(edgesfiltered.all.nocoh$fract_sites_IBD)] <- 1
edgesfiltered.all.nocoh$fract_sites_IBD[is.na(edgesfiltered.all.nocoh$fract_sites_IBD)] <- 0

if(removecontinf){
  print("Removing oversampling of continuous infections")
  edgesfiltered.all.nocoh.nocontinf <- merge(edgesfiltered.all.nocoh, edgesfiltered.nocontinf[,c("ID1", "ID2", "Continf")],
                                             all.x = TRUE)
  edgesfiltered.all.nocoh.nocontinf$Continf[is.na(edgesfiltered.all.nocoh.nocontinf$Continf)] <- FALSE
  edgesfiltered.all.nocoh <- edgesfiltered.all.nocoh.nocontinf[!edgesfiltered.all.nocoh.nocontinf$Continf, -c(ncol(edgesfiltered.all.nocoh.nocontinf))]
  rm(edgesfiltered.all.nocoh.nocontinf)
  nodesfiltered.all.nocoh <- nodesfiltered.all.nocoh[nodesfiltered.all.nocoh$ID%in%c(edgesfiltered.all.nocoh$ID1, 
                                                                                     edgesfiltered.all.nocoh$ID2),]
}


#########VILLAGES##############
minIBDlinks.villages <- 10

nodes.grouped.village <- nodesfiltered.all.nocoh %>%
  group_by(village)%>%
  summarize(barcodes=length(ID), individuals=length(individual))

edges.grouped.village <- edgesfiltered.all.nocoh%>%
  group_by(commonVillage) %>%
  summarize(verticalMax=length(which(sameIndividual)),
            horizontalMax=length(which(!sameIndividual)))

edges.grouped.village.filtered <- edgesfiltered.all.nocoh%>%
  group_by(commonVillage) %>%
  summarize(vertical=length(which(sameIndividual & fract_sites_IBD>=0.5)), 
            horizontal=length(which(!sameIndividual & fract_sites_IBD>=0.5)),
            horizontalIBD=sum(fract_sites_IBD[!sameIndividual & fract_sites_IBD>=0.5]))

edges.grouped.village <- merge(edges.grouped.village.filtered, edges.grouped.village, all.y=T)
edges.grouped.village <- edges.grouped.village[edges.grouped.village$horizontalMax >= minIBDlinks.villages,]
############################

##################COMPOUND################
minIBDlinks.compounds <- 5

nodes.grouped.compound <- nodesfiltered.all.nocoh %>%
  group_by(village, compound)%>%
  summarize(barcodes=length(ID), individuals=length(individual))

edges.grouped.compound <- edgesfiltered.all.nocoh%>%
  group_by(commonVillage, sameVillage, commonCompound, sameCompound) %>%
  summarize(verticalMax=length(which(sameIndividual)),
            horizontalMax=length(which(!sameIndividual)),
            )

edges.grouped.compound.filtered <- edgesfiltered.all.nocoh%>%
  group_by(commonVillage, sameVillage, commonCompound, sameCompound) %>%
  summarize(vertical=length(which(sameIndividual & fract_sites_IBD>=0.5)), 
            horizontal=length(which(!sameIndividual & fract_sites_IBD>=0.5)),
            horizontalIBD=sum(fract_sites_IBD[!sameIndividual & fract_sites_IBD>=0.5]))

edges.grouped.compound <- merge(edges.grouped.compound.filtered, edges.grouped.compound, all.y = TRUE)

edges.grouped.compound <- edges.grouped.compound[edges.grouped.compound$horizontalMax >= minIBDlinks.compounds,]
############################

##################DATE################
minIBDlinks.dates <- 10

nodes.grouped.date <- nodesfiltered.all.nocoh %>%
  group_by(date)%>%
  summarize(barcodes=length(ID), individuals=length(individual))

edges.grouped.date <- edgesfiltered.all.nocoh%>%
  group_by(commonDate, cycleElapsed, sameSeason, commonSeason) %>%
  summarize(verticalMax=length(which(sameIndividual)),
            horizontalMax=length(which(!sameIndividual)),
            datemin = min(date1, date2),
            datemax = max(date1, date2),
            dayElapsed=unique(dayElapsed)
  )

edges.grouped.date.filtered <- edgesfiltered.all.nocoh%>%
  group_by(commonDate, cycleElapsed, sameSeason, commonSeason) %>%
  summarize(vertical=length(which(sameIndividual & fract_sites_IBD>=0.5)), 
            horizontal=length(which(!sameIndividual & fract_sites_IBD>=0.5)),
            horizontalIBD=sum(fract_sites_IBD[!sameIndividual & fract_sites_IBD>=0.5]))

edges.grouped.date <- merge(edges.grouped.date.filtered, edges.grouped.date, all.y = TRUE)

edges.grouped.date <- edges.grouped.date[edges.grouped.date$horizontalMax >= minIBDlinks.dates,]

edges.grouped.date$group <- as.numeric(cut_number(round(edges.grouped.date$dayElapsed/30), 6))
dfgrouplocation <- edges.grouped.date %>%
  group_by(group)%>%
  summarise(mingrouplocation = min(round(dayElapsed/30)),
            maxgrouplocation = max(round(dayElapsed/30)),
            meanIBD = mean(horizontalIBD/horizontalMax))
dfgrouplocation$maxgrouplocation <- c(dfgrouplocation$mingrouplocation[-1], 
                                      max(dfgrouplocation$maxgrouplocation[-1]))
dfgrouplocation$mingrouplocation <- c(min(dfgrouplocation$mingrouplocation), 
                                      dfgrouplocation$maxgrouplocation[-nrow(dfgrouplocation)])
dfgrouplocation$grouplocation <- (dfgrouplocation$maxgrouplocation + dfgrouplocation$mingrouplocation)/2

edges.grouped.date <- merge(edges.grouped.date, dfgrouplocation, by ="group")

morethanonyear <- edges.grouped.date[edges.grouped.date$dayElapsed/30>11.5,]
quantile(morethanonyear$horizontalIBD/morethanonyear$horizontalMax, c(0.05, 0.5, 0.95))
mean(morethanonyear$horizontalIBD/morethanonyear$horizontalMax)

oneseasonapart <- (edges.grouped.date$cycleElapsed <= 1 & 
                     !edges.grouped.date$sameSeason & 
                     abs(year(edges.grouped.date$datemax) - year(edges.grouped.date$datemin)) <= 1) | 
  (edges.grouped.date$cycleElapsed == 0 & edges.grouped.date$sameSeason)

edges.grouped.date.1season <- edges.grouped.date[oneseasonapart,]
edges.grouped.date.not1season <- edges.grouped.date[!oneseasonapart,]
edges.grouped.date.1season%>%
  group_by(commonSeason, cycleElapsed)%>%
  summarise(mean(horizontalIBD/horizontalMax))

drywetwetdry.test <- t.test(horizontalIBD/horizontalMax~paste(commonSeason, cycleElapsed),
                         data=edges.grouped.date.1season[edges.grouped.date.1season$commonSeason=="drywet",])
drywetdry.test <- t.test(horizontalIBD/horizontalMax~paste(commonSeason, cycleElapsed),
                         data=edges.grouped.date.1season[(edges.grouped.date.1season$commonSeason=="drywet" & edges.grouped.date.1season$cycleElapsed==1) | edges.grouped.date.1season$commonSeason=="dry",])
drywetwet.test <- t.test(horizontalIBD/horizontalMax~paste(commonSeason, cycleElapsed),
                         data=edges.grouped.date.1season[(edges.grouped.date.1season$commonSeason=="drywet" & edges.grouped.date.1season$cycleElapsed==1) | edges.grouped.date.1season$commonSeason=="wet",])
drywettests <- data.frame("from" = c("drywet1", "drywet1", "drywet1"))
drywettests$to <- c("dry0", "wet0", "drywet0")
drywettests$tests <- c(drywetwetdry.test$statistic, drywetdry.test$statistic, drywetwet.test$statistic)
drywettests$pvalues <- c(drywetwetdry.test$p.value, drywetdry.test$p.value, drywetwet.test$p.value)

drywet.mean <- edges.grouped.date.1season%>%
  group_by(commonSeason, cycleElapsed)%>%
  summarise(mean(horizontalIBD/horizontalMax))
colnames(drywet.mean) <- c("Season", "Overlap", "Mean")
############################

##################DATE AND COMPOUND################
minIBDlinks.datecompounds <- 5

nodes.grouped.datecompound <- nodesfiltered.all.nocoh %>%
  group_by(date)%>%
  summarize(barcodes=length(ID), individuals=length(individual))

edges.grouped.datecompound <- edgesfiltered.all.nocoh%>%
  group_by(commonDate, cycleElapsed, sameSeason, commonSeason,
           sameVillage, sameCompound) %>%
  summarize(verticalMax=length(which(sameIndividual)),
            horizontalMax=length(which(!sameIndividual)),
            datemin = min(date1, date2),
            datemax = max(date1, date2),
            dayElapsed=unique(dayElapsed)
  )

edges.grouped.datecompound.filtered <- edgesfiltered.all.nocoh%>%
  group_by(commonDate, cycleElapsed, sameSeason, commonSeason,
           sameVillage, sameCompound) %>%
  summarize(vertical=length(which(sameIndividual & fract_sites_IBD>=0.5)), 
            horizontal=length(which(!sameIndividual & fract_sites_IBD>=0.5)),
            horizontalIBD=sum(fract_sites_IBD[!sameIndividual & fract_sites_IBD>=0.5]))

edges.grouped.datecompound <- merge(edges.grouped.datecompound.filtered, edges.grouped.datecompound, all.y = TRUE)

edges.grouped.datecompound <- edges.grouped.datecompound[edges.grouped.datecompound$horizontalMax >= minIBDlinks.datecompounds,]

edges.grouped.datecompound$group <- as.numeric(cut_number(round(edges.grouped.datecompound$dayElapsed/30), 6))
dfgrouplocation.datecompound <- edges.grouped.datecompound %>%
  group_by(group)%>%
  summarise(mingrouplocation = min(round(dayElapsed/30)),
            maxgrouplocation = max(round(dayElapsed/30)),
            meanIBD = mean(horizontalIBD/horizontalMax))
dfgrouplocation.datecompound$maxgrouplocation <- c(dfgrouplocation.datecompound$mingrouplocation[-1], 
                                      max(dfgrouplocation.datecompound$maxgrouplocation[-1]))
dfgrouplocation.datecompound$mingrouplocation <- c(min(dfgrouplocation.datecompound$mingrouplocation), 
                                      dfgrouplocation.datecompound$maxgrouplocation[-nrow(dfgrouplocation.datecompound)])
dfgrouplocation.datecompound$grouplocation <- (dfgrouplocation.datecompound$maxgrouplocation + dfgrouplocation.datecompound$mingrouplocation)/2

edges.grouped.datecompound <- merge(edges.grouped.datecompound, dfgrouplocation.datecompound, by ="group")
############################

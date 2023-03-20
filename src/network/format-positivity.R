#!/usr/bin/env Rscript
library(reshape2)
options(stringsAsFactors = F)
posneg <- read.table("../../../DataEpidemio/Completedata3phases.csv", sep = ";", header = TRUE)
posneg$ParticipantID <- substr(posneg$ParticipantID, 1, 8)

duplIDS <- names(which(table(posneg$ParticipantID)>1))
bestHits <- c()
for (duplID in duplIDS){
  bestHit <- names(which.min(apply(posneg[posneg$ParticipantID==duplID,],1,
                                   FUN=function(x){
                                     length(which(is.na(x)))})))
  bestHits <- c(bestHits, bestHit)
}
posneg <- posneg[c(as.numeric(bestHits), which(!posneg$ParticipantID%in%duplIDS)),]
posneg <- posneg[grepl("^[A-Z]\\d{7}$", posneg$ParticipantID),]
posneg <- posneg[!apply(posneg, 1, FUN= function(x){
  all(is.na(x[-1]))
}),]
length(unique(posneg$ParticipantID))==nrow(posneg)

posneg <- melt(posneg, id.vars = "ParticipantID", 
               variable.name = "date", 
               value.name = "infectivity")
posneg <- posneg[!is.na(posneg$infectivity),]
datemapping <- list("DEC..14"="2014-12-01", "APR..15"="2015-04-01",
                    "JUN..15"="2015-06-01", "NOV..15"="2015-11-01",
                    "DEC..15"="2015-12-01", "MAR..16"="2016-03-01",
                    "MAY..16"="2016-05-01", "EITHERCALL_1607"="2016-07-01",
                    "EITHERCALL_1610"="2016-10-01", "EITHERCALL_1611"="2016-11-01",
                    "CALLATS_1612"="2016-12-01", "DECV.CONSENSUS"="2016-12-15",
                    "JAN.CONSENSUS"="2017-01-01", "FEB.CONSENSUS"="2017-02-01",
                    "MAR.CONSENSUS"="2017-03-01", "APR.CONSENSUS"="2017-04-01",
                    "MAY.CONSENSUS"="2017-05-01")
length(datemapping)==length(unique(toupper(posneg$date)))

posneg$date <- unlist(sapply(posneg$date, FUN=function(x){
  datemapping[[toupper(x)]]
}))
posneg$date <- as.character(as.Date(posneg$date, "%Y-%m-%d"))

posneg$sampleID <- paste(posneg$ParticipantID, format.Date(posneg$date, "%y%m"), sep ="_")

write.csv(posneg[,c("ParticipantID", "sampleID", "date", "infectivity")], 
          "rawdata/pf-test.csv", quote = F, row.names = F)
write.csv(posneg[,c("ParticipantID", "sampleID", "date", "infectivity")], 
          "../../pf-test.csv", quote = F, row.names = F)

###############FORMAT TREATMENT DATA###############
ttm <- read.table("../../../DataEpidemio/Cohorte_true_V1.5_THIS.csv", sep = ";", header = TRUE)
ttm <- ttm[,c("Participant.unique", "Call_1", "Call_2", "Call_3"),]
ttm <- melt(ttm, id.vars = "Participant.unique")
ttm$variable <- as.character(ttm$variable)
ttm$variable[ttm$variable=="Call_1"] <- "2016-07-15"
ttm$variable[ttm$variable=="Call_2"] <- "2016-10-15"
ttm$variable[ttm$variable=="Call_3"] <- "2016-11-15"
ttm$value[ttm$value!=""] <- "Yes"
ttm$value[ttm$value==""] <- NA
colnames(ttm) <- c("ParticipantID", "date", "treatment")
dates <- as.Date(ttm$date, "%Y-%m-%d")

ttm$sampleID <- paste0(ttm$ParticipantID, "_", 
                            substr(year(dates), 3,4),
                            sprintf("%02d", month(dates)))

write.csv(ttm[,c("ParticipantID", "sampleID", "date", "treatment")], 
          "rawdata/pf-treatment.csv", quote = F, row.names = F)
write.csv(ttm[,c("ParticipantID", "sampleID", "date", "treatment")], 
          "../../pf-treatment.csv", quote = F, row.names = F)

####################################

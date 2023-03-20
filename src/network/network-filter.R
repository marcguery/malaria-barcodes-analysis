#################FILTERING#################
filt <- 0.5
datafiltered <- databarcodes[databarcodes$fract_sites_IBD > filt,]

##################################

#############NODES############
nodesfiltered <- datafiltered[!duplicated(datafiltered$ID1),c("ID1", "individual1", "date1", "village1", "compound1")]
colnames(nodesfiltered) <- c("ID2", "individual2", "date2", "village2", "compound2")
nodesfiltered <- rbind(nodesfiltered, datafiltered[!duplicated(datafiltered$sample2),c("ID2", "individual2", "date2", "village2", "compound2")])
nodesfiltered <- nodesfiltered[!duplicated(nodesfiltered$ID2),]

colnames(nodesfiltered) <- c("ID", "individual", "date", "village", "compound")

#################EDGES#################
edgesfiltered <- datafiltered[,c("ID1","ID2", 
                         "fract_sites_IBD", "dayElapsed",
                         "groupElapsed", "sameIndividual", "commonIndividual", 
                         "individual1", "individual2",
                         "sameDate", "commonDate", "date1", "date2",
                         "sameVillage", "commonVillage", "village1", "village2",
                         "sameCompound", "commonCompound", "compound1", "compound2")]
edgesfiltered <- edgesfiltered[with(edgesfiltered, order(date1, date2)),]

combos <- strsplit(paste(edgesfiltered$individual1[!edgesfiltered$sameIndividual], edgesfiltered$individual2[!edgesfiltered$sameIndividual], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edgesfiltered$commonIndividual[!edgesfiltered$sameIndividual] <- combos

combos <- strsplit(paste0(edgesfiltered$village1[!edgesfiltered$sameVillage], edgesfiltered$village2[!edgesfiltered$sameVillage]), split="")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edgesfiltered$commonVillage[!edgesfiltered$sameVillage] <- combos

combos <- strsplit(paste(edgesfiltered$compound1[!edgesfiltered$sameCompound], edgesfiltered$compound2[!edgesfiltered$sameCompound], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edgesfiltered$commonCompound[!edgesfiltered$sameCompound] <- combos

combos <- strsplit(paste(edgesfiltered$date1[!edgesfiltered$sameDate], edgesfiltered$date2[!edgesfiltered$sameDate], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edgesfiltered$commonDate[!edgesfiltered$sameDate] <- combos
##################################

if (!exists("clusterIBDmin")){
  clusterIBDmin <- 0.5
}

edges.filtered <- edgesfiltered[edgesfiltered$fract_sites_IBD >= clusterIBDmin,]
# Get rid of loops and ensure right naming of vertices
edges.topology <- edges.filtered[,c("ID1", "ID2")]
edges.topology <- simplify(graph.data.frame(edges.topology[order(edges.topology[[1]]),],directed = FALSE))

# Find all components
comps <- components(edges.topology)
cluster.df <- data.frame(comps$membership)
cluster.df <- cbind(row.names(cluster.df), cluster.df)
colnames(cluster.df) <- c("ID", "cluster")
nodesfiltered <- merge(nodesfiltered, cluster.df, by = "ID", all.x = T)
nodesfiltered$cluster[is.na(nodesfiltered$cluster)] <- 0

nodesfiltered$cluster <- sprintf(paste0("%0",max(nchar(nodesfiltered$cluster)),"d"), nodesfiltered$cluster)

#################EDGES#################
edgesfiltered <- merge(edgesfiltered, nodesfiltered[,c("ID", "cluster")],
                       by.x = "ID1", by.y = "ID", all.x = TRUE)
edgesfiltered <- merge(edgesfiltered, nodesfiltered[,c("ID", "cluster")],
                       by.x = "ID2", by.y = "ID", all.x = TRUE,
                       suffixes = c("1", "2"))

edgesfiltered$sameCluster <- edgesfiltered$cluster1 == edgesfiltered$cluster2
edgesfiltered$commonCluster <- NA
edgesfiltered$commonCluster[edgesfiltered$sameCluster] <- edgesfiltered$cluster1[edgesfiltered$sameCluster]

combos <- strsplit(paste(edgesfiltered$cluster1[!edgesfiltered$sameCluster], edgesfiltered$cluster2[!edgesfiltered$sameCluster], sep =";"), split=";")
combos <- sapply(combos, sort)
if (length(combos) > 0){
  combos <- paste0(combos[1,], combos[2,])
  edgesfiltered$commonCluster[!edgesfiltered$sameCluster] <- combos
}

##################################

nodesfiltered$date <- as.Date(paste0(nodesfiltered$date, "01"), "%y%m%d")
edgesfiltered$date1 <- as.Date(paste0(edgesfiltered$date1, "01"), "%y%m%d")
edgesfiltered$date2 <- as.Date(paste0(edgesfiltered$date2, "01"), "%y%m%d")

##############SEASONS##############

c_trans <- function(a, b, breaks = b$breaks, format = b$format) {
  a <- as.trans(a)
  b <- as.trans(b)
  
  name <- paste(a$name, b$name, sep = "-")
  
  trans <- function(x) a$trans(b$trans(x))
  inv <- function(x) b$inverse(a$inverse(x))
  
  trans_new(name, trans, inverse = inv, breaks = breaks, format=format)
  
}
rev_date <- c_trans("reverse", "time")
ori_date <- c_trans("identity", "time")

#dry (january to july) wet (beginning august to end december)
wetseason <- c(8,12)
wetseason <- as.character(wetseason)

wetseason[nchar(wetseason)==1] <- paste0("0", wetseason[nchar(wetseason)==1])

years <- sort(as.numeric(unique(format(nodes$date,"%Y"))))
wetseason <- as.Date(as.character(as.Date(c(paste0(years[1], "-01-01"),
                                            unlist(lapply(years, function(x){
                                              paste(x, wetseason, "15", sep="-")
                                            })),
                                            paste0(years[length(years)], "-12-31")), "%Y-%m-%d")), "%Y-%m-%d")

seasons <- data.frame(wetseason[-length(wetseason)], wetseason[-1])
seasons$type <- rep(c("dry", "wet"), length.out=nrow(seasons))
seasons$cycle <- rep(c(-1,-1), length.out=nrow(seasons))
seasons$cycle <- ifelse(rep(seasons$type[1]=="wet", length.out=nrow(seasons)),
                        seasons$cycle+floor((0:(nrow(seasons)-1))/2)+1,
                        seasons$cycle+floor((1:nrow(seasons))/2)+1)
colnames(seasons) <- c("date1", "date2", "type", "cycle")
mindate <- min(nodes$date)
maxdate <- max(nodes$date)
#Time zone warning but it is ok
seasons <- seasons[seasons$date2>mindate & seasons$date1 < maxdate,]

seasons.combination <- expand.grid(paste(seasons$date1, seasons$date2),
                                   paste(seasons$date1, seasons$date2), stringsAsFactors = F)
seasons.combination$date1start <- as.POSIXct(as.character(as.Date(sub("\\s\\S+", "", seasons.combination[,1]), "%Y-%m-%d")), "%Y-%m-%d")
seasons.combination$date1end <- as.POSIXct(as.character(as.Date(sub("\\S+\\s", "", seasons.combination[,1]), "%Y-%m-%d")), "%Y-%m-%d")
seasons.combination$date2start <- as.POSIXct(as.character(as.Date(sub("\\s\\S+", "", seasons.combination[,2]), "%Y-%m-%d")), "%Y-%m-%d")
seasons.combination$date2end <- as.POSIXct(as.character(as.Date(sub("\\S+\\s", "", seasons.combination[,2]), "%Y-%m-%d")), "%Y-%m-%d")
seasons.combination <- seasons.combination[,-c(1,2)]

seasons.combination <- merge(seasons.combination, seasons,
                             by.x=c("date1start", "date1end"), by.y=c("date1", "date2"))
colnames(seasons.combination)[c((ncol(seasons.combination)-1):ncol(seasons.combination))] <- c("season1", "cycle1")
seasons.combination <- merge(seasons.combination, seasons, 
                             by.x=c("date2start", "date2end"), by.y=c("date1", "date2"))
colnames(seasons.combination)[c((ncol(seasons.combination)-1):ncol(seasons.combination))] <- c("season2", "cycle2")
############################

nodesfiltered$season <- apply(nodesfiltered, 1, 
                              FUN = function(x){
                                seasons$type[seasons$date1 <= x [3] & seasons$date2 >= x[3]]})
nodesfiltered$cycle <- apply(nodesfiltered, 1, 
                             FUN = function(x){
                               seasons$cycle[seasons$date1 <= x [3] & seasons$date2 >= x[3]]})




  edgesfiltered <- merge(edgesfiltered, nodesfiltered[,c("ID", "season", "cycle")], all.x = T,
                         by.x = "ID1", by.y = "ID")
  edgesfiltered <- merge(edgesfiltered, nodesfiltered[,c("ID", "season", "cycle")], all.x = T,
                         by.x = "ID2", by.y = "ID", suffixes = c("1", "2"))
  


##Seasons
edgesfiltered$sameSeason <- edgesfiltered$season1==edgesfiltered$season2
edgesfiltered$commonSeason <- NA
edgesfiltered$commonSeason[edgesfiltered$sameSeason] <- edgesfiltered$season1[edgesfiltered$sameSeason]

combos <- strsplit(paste(edgesfiltered$season1[!edgesfiltered$sameSeason], edgesfiltered$season2[!edgesfiltered$sameSeason], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edgesfiltered$commonSeason[!edgesfiltered$sameSeason] <- combos

##Cycles
edgesfiltered$sameCycle <- edgesfiltered$cycle1==edgesfiltered$cycle2
edgesfiltered$commonCycle <- NA
edgesfiltered$commonCycle[edgesfiltered$sameCycle] <- edgesfiltered$cycle1[edgesfiltered$sameCycle]

combos <- strsplit(paste(edgesfiltered$cycle1[!edgesfiltered$sameCycle], edgesfiltered$cycle2[!edgesfiltered$sameCycle], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edgesfiltered$commonCycle[!edgesfiltered$sameCycle] <- combos

edgesfiltered$cycleElapsed <- abs(edgesfiltered$cycle2 - edgesfiltered$cycle1)

###########################SAVE DATA###########################

write.table(datafiltered, 
            file = paste0("out/filtered-barcodes-", 
                          sub("\\.", "", as.character(filt)), 
                          "-all.hmm_fract.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(nodesfiltered,
            file = "out/filtered-nodes.tsv",
            quote = F, row.names = F)
write.table(edgesfiltered,
            file = "out/filtered-edges.tsv",
            quote = F, row.names = F)
######################################################




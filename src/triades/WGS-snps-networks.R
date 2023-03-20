###########################HMM IBD FILE###########################
hmmfile <- "../hmmibd/out/IBD-WGS-snps.hmm_fract.txt.zip"
databarcodes <- read_delim(file = hmmfile, delim = "\t")
barcodes <- read.csv("../read/out/WGS-19262snps_barcodes.csv")
samplestoremove <- read.table("../read/rawdata/samples-to-remove.txt", sep ="\t", header = TRUE)
######################################################

################WGS snps MAPPING#################
WGSmapping <- read.csv("../read/out/WGS-19262snps_mapping.csv")
colnames(WGSmapping) <- c(colnames(WGSmapping)[-ncol(WGSmapping)], "ID")

barcodes <- merge(barcodes, WGSmapping[,c(1,3)], 
                     by.x = "Internal.Sample.ID", by.y = "sample", all.x = T)
barcodes$Internal.Sample.ID <- make.names(barcodes$Internal.Sample.ID)
barcodes$ID[barcodes$Internal.Sample.ID%in%c("REF.WGS", "ALT.WGS")] <- barcodes$Internal.Sample.ID[barcodes$Internal.Sample.ID%in%c("REF.WGS", "ALT.WGS")]
barcodes <- barcodes[!is.na(barcodes$ID),]
################################

###########################ATTRIBUTING IDs###########################
databarcodes <- merge(databarcodes, barcodes[,c("Internal.Sample.ID", "ID")], 
                      by.x = "sample1", by.y = "Internal.Sample.ID")
colnames(databarcodes)[ncol(databarcodes)] <- "ID1"
databarcodes <- merge(databarcodes, barcodes[,c("Internal.Sample.ID", "ID")], 
                      by.x = "sample2", by.y = "Internal.Sample.ID")
colnames(databarcodes)[ncol(databarcodes)] <- "ID2"

all(table(c(databarcodes$sample1, databarcodes$sample2))==201)
all(table(c(databarcodes$ID1, databarcodes$ID2))==201)

minsites <- 100
databarcodes$fract_sites_IBD[databarcodes$N_informative_sites<minsites] <- -1
databarcodes$seq_shared_best_traj[databarcodes$N_informative_sites<minsites] <- -1
######################################################

###########################FILTERING OPTIONS###########################
#Removal of 3D7 cluster: removing nodes having IBD >= 0.9 with REF
#This does not remove (J0171806) which is indirectly linked to REF with IBD > 0.9

databarcodes.ref <- databarcodes[databarcodes$fract_sites_IBD>=0.9 & (grepl("REF.WGS", databarcodes$sample1) | grepl("REF.WGS", databarcodes$sample2)),]
samplesref <- unique(c(databarcodes.ref$sample1, databarcodes.ref$sample2, "REF.WGS"))
databarcodes <- databarcodes[!databarcodes$sample1%in%samplesref & !databarcodes$sample2%in%samplesref,]
databarcodes.alt <- databarcodes[databarcodes$fract_sites_IBD>=0.9 & (grepl("ALT.WGS", databarcodes$sample1) | grepl("ALT.WGS", databarcodes$sample2)),]
samplesalt <- unique(c(databarcodes.alt$sample1, databarcodes.alt$sample2, "ALT.WGS"))
databarcodes <- databarcodes[!databarcodes$sample1%in%samplesalt & !databarcodes$sample2%in%samplesalt,]
databarcodes <- databarcodes[!databarcodes$sample1%in%samplestoremove$ID & !databarcodes$sample2%in%samplestoremove$ID,]
######################################################

###########################SPACE/TIME INFORMATION###########################
##DATES
databarcodes$date1 <- str_extract(databarcodes$ID1, pattern = "[0-9]{4}$")
databarcodes$date2 <- str_extract(databarcodes$ID2, pattern = "[0-9]{4}$")
date1 <- databarcodes$date1
date2 <- databarcodes$date2

date1[!is.na(date1)] <- paste(date1[!is.na(date1)], "01", sep="")
date2[!is.na(date2)] <- paste(date2[!is.na(date2)], "01", sep="")
date1.date <- as.Date(date1, "%y%m%d")
date2.date <- as.Date(date2, "%y%m%d")
diffdate <- round(abs(difftime(date1.date,date2.date, units="days")), digits = 2)
dateorder <- abs(as.numeric(as.factor(date1))-as.numeric(as.factor(date2)))
databarcodes$dayElapsed <- diffdate
databarcodes$groupElapsed <- dateorder
databarcodes$sameDate <- date1==date2
databarcodes$commonDate <- NA
databarcodes$commonDate[databarcodes$sameDate] <- databarcodes$date1[databarcodes$sameDate]

##LOCATIONS
databarcodes$individual1 <- sub(pattern = "_\\S+", replacement = "", x = databarcodes$ID1)
databarcodes$individual2 <- sub(pattern = "_\\S+", replacement = "", x = databarcodes$ID2)
databarcodes$sameIndividual <- databarcodes$individual1==databarcodes$individual2
databarcodes$commonIndividual <- NA
databarcodes$commonIndividual[databarcodes$sameIndividual] <- databarcodes$individual1[databarcodes$sameIndividual]


databarcodes$village1 <- substr(databarcodes$ID1, 1, 1)
databarcodes$village2 <- substr(databarcodes$ID2, 1, 1)
databarcodes$sameVillage <- databarcodes$village1==databarcodes$village2
databarcodes$commonVillage <- NA
databarcodes$commonVillage[databarcodes$sameVillage] <- databarcodes$village1[databarcodes$sameVillage]

databarcodes$compound1 <- substr(databarcodes$ID1, 1, 4)
databarcodes$compound2 <- substr(databarcodes$ID2, 1, 4)
databarcodes$sameCompound <- databarcodes$compound1==databarcodes$compound2
databarcodes$commonCompound <- NA
databarcodes$commonCompound[databarcodes$sameCompound] <- databarcodes$compound1[databarcodes$sameCompound]

#Number of samples
length(unique(c(databarcodes$ID1, databarcodes$ID2)))
#Number of individuals
length(unique(c(databarcodes$individual1, databarcodes$individual2)))

databarcodes.all <- databarcodes
######################################################

#############NODES############
nodes <- databarcodes[!duplicated(databarcodes$ID1),c("ID1", "individual1", "date1", "village1", "compound1")]
colnames(nodes) <- c("ID2", "individual2", "date2", "village2", "compound2")
nodes <- rbind(nodes, databarcodes[!duplicated(databarcodes$sample2),c("ID2", "individual2", "date2", "village2", "compound2")])
nodes <- nodes[!duplicated(nodes$ID2),]

colnames(nodes) <- c("ID", "individual", "date", "village", "compound")
##########################
#############FIND PARENTS#############
edges.filtered <- databarcodes[databarcodes$seq_shared_best_traj>=parentalmin & databarcodes$seq_shared_best_traj<=parentalmax,]
# Get rid of loops and ensure right naming of vertices
edges.topology <- edges.filtered[,c("ID1", "ID2")]
edges.topology <- simplify(graph.data.frame(edges.topology[order(edges.topology[[1]]),],directed = FALSE))

# Find all components
parent <- components(edges.topology)
parent.df <- data.frame(parent$membership)
parent.df <- cbind(row.names(parent.df), parent.df)
colnames(parent.df) <- c("ID", "parent")
nodes <- merge(nodes, parent.df, all.x = T)
nodes$parent[which(is.na(nodes$parent))] <- c(1:length(which(is.na(nodes$parent))))+max(nodes$parent, na.rm = TRUE)

nodes$parent <- sprintf(paste0("%0",max(nchar(nodes$parent)),"d"), nodes$parent)
##########################

#################FIND IDENTICAL SAMPLES#################
edges.filtered <- databarcodes[databarcodes$seq_shared_best_traj>=identical,]
# Get rid of loops and ensure right naming of vertices
edges.topology <- edges.filtered[,c("ID1", "ID2")]
edges.topology <- simplify(graph.data.frame(edges.topology[order(edges.topology[[1]]),],directed = FALSE))

# Find all components
strain <- components(edges.topology)
strain.df <- data.frame(strain$membership)
strain.df <- cbind(row.names(strain.df), strain.df)
colnames(strain.df) <- c("ID", "strain")
nodes <- merge(nodes, strain.df, all.x = T)
nodes$strain[which(is.na(nodes$strain))] <- c(1:length(which(is.na(nodes$strain))))+max(nodes$strain, na.rm = TRUE)


nodes$strain <- sprintf(paste0("%0",max(nchar(nodes$strain)),"d"), nodes$strain)
##########################
#################EDGES#################
edges <- databarcodes[,c("ID1","ID2", 
                         "fract_sites_IBD", "seq_shared_best_traj", "dayElapsed",
                         "groupElapsed", "sameIndividual", "commonIndividual", 
                         "individual1", "individual2",
                         "sameDate", "commonDate", "date1", "date2",
                         "sameVillage", "commonVillage", "village1", "village2",
                         "sameCompound", "commonCompound", "compound1", "compound2")]
edges <- edges[with(edges, order(date1, date2)),]

combos <- strsplit(paste(edges$individual1[!edges$sameIndividual], edges$individual2[!edges$sameIndividual], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonIndividual[!edges$sameIndividual] <- combos

combos <- strsplit(paste0(edges$village1[!edges$sameVillage], edges$village2[!edges$sameVillage]), split="")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonVillage[!edges$sameVillage] <- combos

combos <- strsplit(paste(edges$compound1[!edges$sameCompound], edges$compound2[!edges$sameCompound], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonCompound[!edges$sameCompound] <- combos

combos <- strsplit(paste(edges$date1[!edges$sameDate], edges$date2[!edges$sameDate], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonDate[!edges$sameDate] <- combos

edges <- merge(edges, nodes[,c("ID", "parent", "strain")],
               by.x = "ID1", by.y = "ID", all.x = TRUE)
edges <- merge(edges, nodes[,c("ID", "parent", "strain")],
               by.x = "ID2", by.y = "ID", all.x = TRUE,
               suffixes = c("1", "2"))

edges$sameParent <- edges$parent1 == edges$parent2
edges$commonParent <- NA
edges$commonParent[edges$sameParent] <- edges$parent1[edges$sameParent]

combos <- strsplit(paste(edges$parent1[!edges$sameParent], edges$parent2[!edges$sameParent], sep =";"), split=";")
combos <- sapply(combos, sort)
if (length(combos) > 0){
  combos <- paste0(combos[1,], combos[2,])
  edges$commonParent[!edges$sameParent] <- combos
}

edges$sameStrain <- edges$strain1 == edges$strain2
edges$commonStrain <- NA
edges$commonStrain[edges$sameStrain] <- edges$strain1[edges$sameStrain]

combos <- strsplit(paste(edges$strain1[!edges$sameStrain], edges$strain2[!edges$sameStrain], sep =";"), split=";")
combos <- sapply(combos, sort)
if (length(combos) > 0){
  combos <- paste0(combos[1,], combos[2,])
  edges$commonStrain[!edges$sameStrain] <- combos
}

##################################

nodes$date <- as.Date(paste0(nodes$date, "01"), "%y%m%d")
edges$date1 <- as.Date(paste0(edges$date1, "01"), "%y%m%d")
edges$date2 <- as.Date(paste0(edges$date2, "01"), "%y%m%d")

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

nodes$season <- apply(nodes, 1, 
                      FUN = function(x){
                        seasons$type[seasons$date1 <= x [3] & seasons$date2 >= x[3]]})
nodes$cycle <- apply(nodes, 1, 
                     FUN = function(x){
                       seasons$cycle[seasons$date1 <= x [3] & seasons$date2 >= x[3]]})


edges <- merge(edges, nodes[,c("ID", "season", "cycle")], all.x = T,
               by.x = "ID1", by.y = "ID")
edges <- merge(edges, nodes[,c("ID", "season", "cycle")], all.x = T,
               by.x = "ID2", by.y = "ID", suffixes = c("1", "2"))


##Seasons
edges$sameSeason <- edges$season1==edges$season2
edges$commonSeason <- NA
edges$commonSeason[edges$sameSeason] <- edges$season1[edges$sameSeason]

combos <- strsplit(paste(edges$season1[!edges$sameSeason], edges$season2[!edges$sameSeason], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonSeason[!edges$sameSeason] <- combos

##Cycles
edges$sameCycle <- edges$cycle1==edges$cycle2
edges$commonCycle <- NA
edges$commonCycle[edges$sameCycle] <- edges$cycle1[edges$sameCycle]

combos <- strsplit(paste(edges$cycle1[!edges$sameCycle], edges$cycle2[!edges$sameCycle], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonCycle[!edges$sameCycle] <- combos

edges$cycleElapsed <- abs(edges$cycle2 - edges$cycle1)

###########################STATS###########################
mostsimilars <- edges[edges$fract_sites_IBD>0.99,]
mostsimilars.id <- unique(c(mostsimilars$ID1, mostsimilars$ID2))
mostsimilars <- expand.grid(mostsimilars.id, mostsimilars.id)
colnames(mostsimilars) <- c("ID1", "ID2")
mostsimilars <- merge(mostsimilars, edges)
######################################################

###########################SAVE DATA###########################
write.table(databarcodes, 
            file = paste0("out/WGS-snps-all.hmm_fract.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(nodes,
            file = "out/nodes.tsv",
            quote = F, row.names = F)
write.table(edges,
            file = "out/edges.tsv",
            quote = F, row.names = F)

write.table(mostsimilars,
            file = "out/similar-samples.tsv",
            quote = F, row.names = F)
##################################

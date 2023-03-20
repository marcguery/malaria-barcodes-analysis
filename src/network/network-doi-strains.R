###########################SCRIPT###########################
nodesfilteredclusters <- data.frame(ID = NA, individual = NA, cluster = NA, IBD = NA)
nodesfilteredlineages <- data.frame(individual = NA, cluster = NA, days = NA, 
                            start = NA, end = NA, IBD = NA)
edgesfilteredsubsists <- data.frame(commonCluster = NA, days = NA,
                            start = NA, end = NA, IBD = NA)
nodesfilteredsubsists <- data.frame(ID = NA, individual = NA, date = NA,
                            village = NA, compound = NA, cluster= NA, 
                            season = NA, cycle = NA, IBD = NA)
#All these cutoffs will be tested
clusterIBDmins <- seq(1,0.5,-0.025)
for (clusterIBDmin in clusterIBDmins){
  source("network-filter.R")
  source("duration-of-infection.R")
  nodesfilteredcluster <- nodesfiltered[,c("ID", "individual", "cluster")]
  nodesfilteredcluster$IBD <- clusterIBDmin
  nodesfilteredclusters <- rbind(nodesfilteredclusters, nodesfilteredcluster)
  nodesfilteredlineage <- nodesfiltered.lasting
  nodesfilteredlineage$IBD <- clusterIBDmin
  nodesfilteredlineages <- rbind(nodesfilteredlineages, nodesfilteredlineage)
  edgesfilteredsubsist <- edgesfiltered.subsist
  edgesfilteredsubsist$IBD <- clusterIBDmin
  edgesfilteredsubsists <- rbind(edgesfilteredsubsists, edgesfilteredsubsist)
  nodesfilteredsubsist <- nodesfiltered.subsist
  nodesfilteredsubsist$IBD <- clusterIBDmin
  nodesfilteredsubsists <- rbind(nodesfilteredsubsists, nodesfilteredsubsist)
}


nodesfilteredclusters <- nodesfilteredclusters[!is.na(nodesfilteredclusters$ID),]
nodesfilteredlineages <- nodesfilteredlineages[!is.na(nodesfilteredlineages$individual),]
edgesfilteredsubsists <- edgesfilteredsubsists[!is.na(edgesfilteredsubsists$commonCluster),]
nodesfilteredsubsists <- nodesfilteredsubsists[!is.na(nodesfilteredsubsists$ID),]
###########################
###########################2. NUMBER OF STRAINS###########################

#Epidemiologic data
posneg <- read.csv("rawdata/pf-test.csv")
posneg$date <- as.Date(posneg$date)

ttm <- read.csv("rawdata/pf-treatment.csv")
ttm$date <- as.Date(ttm$date)

posnegttm <- merge(posneg, ttm, by = c("ParticipantID", "sampleID"),
                   suffixes = c("Posneg", "Ttm"), all = T)

#We chose the cutoff maximizing N matches and minimizing N mismacthes
nodesfilteredclusters.summary <- nodesfilteredclusters %>%
  group_by(IBD, cluster) %>%
  summarise(individuals = length(unique(ID)))

nodesfilteredclusters.summary <- nodesfilteredclusters.summary[with(nodesfilteredclusters.summary, order(IBD, individuals, cluster)),]
nodesfilteredclusters.summary$groups <- factor(paste0(nodesfilteredclusters.summary$IBD, nodesfilteredclusters.summary$individuals, nodesfilteredclusters.summary$cluster),
                                       levels = paste0(nodesfilteredclusters.summary$IBD, nodesfilteredclusters.summary$individuals, nodesfilteredclusters.summary$cluster))
gg <- ggplot(nodesfilteredclusters.summary)+
  geom_bar(aes(x = IBD, y = individuals, group = groups,
               fill = cluster != "00"), 
           color = "black", stat = 'identity')
gg

ggsave("out/minibdforclusters.png", width = 8, height = 6)

######################################################

###########################3.DURATION OF INFECTION###########################

#We chose the cutoff maximizing N matches and minimizing N mismacthes
nodesfilteredlineages$start <- as.POSIXct(as.Date(nodesfilteredlineages$start, origin = "1970-01-01"))
nodesfilteredlineages$end <- as.POSIXct(as.Date(nodesfilteredlineages$end, origin = "1970-01-01"))

nodesfilteredlineages <- nodesfilteredlineages[with(nodesfilteredlineages, order(IBD, -days, start)),]
nodesfilteredlineages$groups <- factor(paste0(nodesfilteredlineages$IBD,
                                      nodesfilteredlineages$individual,
                                      nodesfilteredlineages$cluster),
                               levels = paste0(nodesfilteredlineages$IBD,
                                               nodesfilteredlineages$individual,
                                               nodesfilteredlineages$cluster))
posnegttmlineages <- merge(posnegttm, nodesfilteredlineages, 
                              by.x = "ParticipantID", by.y = "individual")

ibdcutoff <- 0.5
posnegttmlineages.subset <- posnegttmlineages[posnegttmlineages$IBD == ibdcutoff,]
posnegttmlineages.subset1 <- posnegttmlineages.subset[!is.na(posnegttmlineages.subset$datePosneg) & (posnegttmlineages.subset$datePosneg + 180 > posnegttmlineages.subset$start) &
                            (posnegttmlineages.subset$datePosneg - 180 < posnegttmlineages.subset$end),]
posnegttmlineages.subset2 <- posnegttmlineages.subset[is.na(posnegttmlineages.subset$datePosneg) & (posnegttmlineages.subset$dateTtm + 180 > posnegttmlineages.subset$start) &
                                                       (posnegttmlineages.subset$dateTtm - 180 < posnegttmlineages.subset$end),]
posnegttmlineages.subset <- rbind(posnegttmlineages.subset1, posnegttmlineages.subset2)

nodesfilteredlineages.subset <- nodesfilteredlineages[nodesfilteredlineages$IBD == ibdcutoff,]
#nodesfilteredlineages.subset$individual[nodesfilteredlineages.subset$individual=="J0282916"] <- "P0060623"

nodesfilteredlineages.subset <- nodesfilteredlineages.subset[with(nodesfilteredlineages.subset, 
                                                                  order(individual, start, cluster)),]


nodesfilteredlineages.subset <- merge(nodesfilteredlineages.subset,nodesfilteredlineages.subset%>%
  group_by(individual)%>%
    arrange(individual, start)%>%
    filter(length(cluster) > 1) %>%
  summarise(groups = groups, 
            rank = as.numeric(factor(cluster, levels = cluster)),
            relpos = seq(0.75,1.25,length.out = length(cluster))),
  all = T)
nodesfilteredlineages.subset[is.na(nodesfilteredlineages.subset$rank),c("rank", "relpos")] <- 1
nodesfilteredlineages.subset <- nodesfilteredlineages.subset[with(nodesfilteredlineages.subset, 
                                                                  order(-days, start)),]
nodesfilteredlineages.subset$individual <- factor(nodesfilteredlineages.subset$individual, 
                                                     levels = unique(nodesfilteredlineages.subset$individual))

seasons.narrow <- seasons
seasons.narrow$date1[seasons.narrow$date1==min(seasons.narrow$date1)] <- "2014-11-15"
seasons.narrow$date2[seasons.narrow$date2==max(seasons.narrow$date2)] <- "2017-05-30"

posnegttmlineages.subset <- subset(posnegttmlineages.subset, !is.na(infectivity))
#posnegttmlineages.subset$ParticipantID[posnegttmlineages.subset$ParticipantID=="J0282916"] <- "P0060623"
posnegttmlineages.subset <- posnegttmlineages.subset[!duplicated(posnegttmlineages.subset[,c("ParticipantID", "datePosneg")]),]

axisshift <- length(unique(nodesfilteredlineages.subset$individual))/20

gg <- ggplot(nodesfilteredlineages.subset)+
  geom_rect(data=seasons.narrow, aes(ymin=-Inf, ymax=Inf, 
                              xmin=as.POSIXct(date1), xmax=as.POSIXct(date2), 
                              fill=toupper(type)), alpha=0.1)+
  scale_fill_manual("Season",values=c("gold2", "darkblue"))+
  geom_linerange(data = posnegttmlineages.subset,
             aes(ymin = as.numeric(factor(ParticipantID, 
                                        levels = levels(nodesfilteredlineages.subset$individual)))+0.325,
                 ymax = as.numeric(factor(ParticipantID, 
                            levels = levels(nodesfilteredlineages.subset$individual)))-0.325,
                 x = as.POSIXct(datePosneg),
                 color = infectivity), 
             alpha = 0.6, size = 0.5)+
  scale_color_manual("Infectivity",
                     breaks = c("NEG", "POS"),
                     labels = c("Pf-", "Pf+"),
                     values = c("blue", "red"),)+
  geom_segment(aes(y=as.numeric(factor(individual))-0.5, 
                   yend=as.numeric(factor(individual))-0.5,
                   x=as.POSIXct("2014-11-15"), 
                   xend=as.POSIXct("2017-05-31")), 
               size = 0.25, color = "grey80")+
  geom_segment(data=data.frame(1),
               y=max(as.numeric(factor(nodesfilteredlineages.subset$individual)))+0.5, 
                   yend=max(as.numeric(factor(nodesfilteredlineages.subset$individual)))+0.5,
                   x=as.POSIXct("2014-11-15"), 
                   xend=as.POSIXct("2017-05-31"), 
               size = 0.25, color = "grey80")+
  geom_segment(aes(y=as.numeric(factor(individual))-(1-relpos), 
                   yend=as.numeric(factor(individual))-(1-relpos),
                   x=start, xend=end, group = rank, linetype = factor(rank)), size = 0.25,
               show.legend = F)+
  
  geom_point(data = posnegttmlineages.subset[!is.na(posnegttmlineages.subset$treatment),],
                 aes(y = as.numeric(factor(ParticipantID, 
                                              levels = levels(nodesfilteredlineages.subset$individual))),
                     x = as.POSIXct(dateTtm)), color = "green4",
             pch = 8, size = 1)+
  geom_point(aes(y = as.numeric(factor(individual))-(1-relpos), 
                 x = start, 
                 group = rank), size = 0.5,
             show.legend = F)+
  geom_point(aes(y = as.numeric(factor(individual))-(1-relpos), 
                 x = end, 
                 group = rank), size = 0.5,
             show.legend = F)+
  scale_x_continuous(trans = ori_date,
                     breaks=as.POSIXct(unique(c(nodesfiltered$date))), 
                     labels = month(unique(c(as.Date(nodesfiltered$date))), 
                                    label = T, abbr = T),
                     expand = c(0,0))+
  scale_y_continuous(breaks = c(1:length(unique(nodesfilteredlineages.subset$individual))),
                     labels = levels(nodesfilteredlineages.subset$individual),
                     expand = c(0.017,0))+
  ylab("")+
  xlab("")+
  annotate("text",
           label = "2014",
           x=as.POSIXct("2014-12-7"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2014-11-18"),
           xend=as.POSIXct("2014-12-30"),
           y=-axisshift+0.75,
           yend=-axisshift+0.75,
           color ="grey60") +
  annotate("text",
           label = "2015",
           x=as.POSIXct("2015-06-22"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2015-01-03"),
           xend=as.POSIXct("2015-12-30"),
           y=-axisshift+0.75,
           yend=-axisshift+0.75,
           color ="grey60") +
  annotate("text",
           label = "2016",
           x=as.POSIXct("2016-06-22"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2016-01-03"),
           xend=as.POSIXct("2016-12-30"),
           y=-axisshift+0.75,
           yend=-axisshift+0.75,
           color ="grey60") +
  annotate("text",
           label = "2017",
           x=as.POSIXct("2017-03-22"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2017-01-03"),
           xend=as.POSIXct("2017-05-30"),
           y=-axisshift+0.75,
           yend=-axisshift+0.75,
           color ="grey60") +
  coord_cartesian(ylim=c(1,length(unique(nodesfilteredlineages.subset$individual))), 
                  clip="off")+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        legend.title= element_text(size = 10),
        legend.margin = margin(0,0,0,0),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(1,1,2,1), "lines"),
        panel.background = element_blank(),
        text = element_text(size = 14))
  
gg
ggsave(paste0("out/doi-", sub("\\.", "", as.character(ibdcutoff), perl = F), ".png"), width = 13, height = 8)


nodes.unique <- data.frame(unique(nodes[,c("individual")]))
colnames(nodes.unique) <- "individual"
nodesfilteredlineages.volunteer <- merge(nodesfilteredlineages.subset, nodes.unique, all = T)

studyParticipants <- read.csv("../../study-participants.csv", header = T)
nodesfilteredlineages.volunteer <- merge(nodesfilteredlineages.volunteer,
                                                studyParticipants,
                                                by.x = "individual", by.y = "SubjectID", all.x = T)

nodesfilteredlineages.volunteer <- nodesfilteredlineages.volunteer[!is.na(nodesfilteredlineages.volunteer$DOB),]

gg <- ggplot(nodesfilteredlineages.volunteer)+
  geom_bar(aes(x = round(DOB/10)*10, 
               fill = factor(paste0(Gender, is.na(days)),
               levels = c("MTRUE", "MFALSE", "FTRUE", "FFALSE"))),
             position = "stack", alpha = 0.8)+
  scale_x_continuous(n.breaks = 10, expand = c(0,0),
                     labels = function(x){paste0(x,"-",x+9)})+
  scale_y_continuous(n.breaks = 10, expand = c(0,0), limits = c(0,95))+
  scale_fill_manual("",
                    breaks = c("FFALSE", "FTRUE", "MFALSE", "MTRUE"),
                     values = c("blue", "lightblue", "red", "pink"),
                     labels = c("Female (DOI > 1 month)","Female (DOI unknown)",
                                "Male (DOI > 1 month)","Male (DOI unknown)"))+
  xlab("Year of birth")+
  ylab("Individuals")+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        text = element_text(size = 18))
  
gg
ggsave("out/doi-avail-yob.png", width = 12, height = 6, dpi = 150)

gg <- ggplot(nodesfilteredlineages.volunteer[!is.na(nodesfilteredlineages.volunteer$DOB) & !is.na(nodesfilteredlineages.volunteer$days),])+
  geom_vline(xintercept = c(1965+seq(0,50,10)),
             linetype = 2)+
  geom_point(aes(x = round(DOB/10)*10+(as.numeric(factor(Gender))-1.5)*5,
                 y = round(days/30), 
                 fill = Gender),
             alpha = 0.8, pch = 21, size = 2, 
             position = position_jitter(height = 0, seed = 123))+
  scale_fill_manual("",
                    breaks = c("F", "M"),
                     values = c("blue", "red"),
                     labels = c("Female", "Male"))+
  xlab("Year of birth")+
  ylab("Duration of infection (months)")+
  scale_x_continuous(expand = c(0,0), limits = c(1964.5, 2015.5),
                     labels = function(x){paste0(x,"-",x+9)})+
  scale_y_continuous(n.breaks = 15, expand = c(0,0), limits = c(0,18))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        text = element_text(size = 18))
gg
ggsave("out/doi-yob.png", width = 8, height = 6, dpi = 150)

nodesfilteredlineages.volunteer$startage <- time_length(difftime(as.Date(nodesfilteredlineages.volunteer$start), as.Date(paste0(nodesfilteredlineages.volunteer$DOB, "-01-01"))), "years")
gg <- ggplot(nodesfilteredlineages.volunteer[!is.na(nodesfilteredlineages.volunteer$DOB) & !is.na(nodesfilteredlineages.volunteer$days),])+
  geom_vline(xintercept = c(5,15,25,35),
             linetype = 2)+
  geom_point(aes(x = pmin(round(startage/10)*10,30)+5*(as.numeric(factor(Gender))-1.5),
                 y = round(days/30), 
                 fill = Gender),
             alpha = 0.8, pch = 21, size = 2, 
             position = position_jitter(height = 0, seed = 123))+
  scale_fill_manual("",
                    breaks = c("F", "M"),
                    values = c("blue", "red"),
                    labels = c("Female", "Male"))+
  xlab("Age")+
  ylab("Duration of infection (months)")+
  scale_x_continuous(breaks = c(10, 20, 30),
                     labels = c("5-15", "15-25", "25+"),
                     expand = c(0,0), limits = c(4.75,35.25))+
  scale_y_continuous(n.breaks = 15, expand = c(0,0), limits = c(0,18))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        text = element_text(size = 18))
gg
ggsave("out/doi-age.png", width = 8, height = 6, dpi = 150)

gg <- ggplot(nodesfilteredlineages.volunteer)+
  geom_bar(aes(x = round(DOB/10)*10, 
               fill = factor(floor(days/60)),
               group = factor(-floor(days/60),
                              levels = c(NA, unique(-floor(days[!is.na(days)]/60))), exclude = NULL)),
           position = "fill")+
  scale_fill_brewer("Duration of infection",
                    type = "seq", palette = "Greens", 
                    na.value = NA,
                    breaks = c("0","1","2","3","5","8"),
                    labels = c("1-2 months", "2-4 months",
                               "4-6 months", "6-10 months",
                               "10-16 months", "> 16 months"))+
  scale_x_continuous(n.breaks = 5, expand = c(0,0), limits = c(1965,2015),
                     labels = function(x){paste0(x,"-",x+9)})+
  scale_y_continuous(n.breaks = 8, expand = c(0,0), limits = c(0,0.31))+
  xlab("Year of birth")+
  ylab("Proportion over positive individuals")+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        text = element_text(size = 18))
gg
ggsave("out/doi-yob-prop.png", width = 10, height = 6, dpi = 150)


nodesfilteredsubsists$date <- as.POSIXct(as.Date(nodesfilteredsubsists$date, origin = "1970-01-01"))
edgesfilteredsubsists$start <- as.POSIXct(as.Date(edgesfilteredsubsists$start, origin = "1970-01-01"))
edgesfilteredsubsists$end <- as.POSIXct(as.Date(edgesfilteredsubsists$end, origin = "1970-01-01"))

edgesfilteredsubsists <- edgesfilteredsubsists[with(edgesfilteredsubsists, order(IBD, days)),]
edgesfilteredsubsists$groups <- factor(paste0(edgesfilteredsubsists$IBD,
                                              edgesfilteredsubsists$commonCluster),
                                       levels = paste0(edgesfilteredsubsists$IBD,
                                                       edgesfilteredsubsists$commonCluster))

nodesfilteredsubsists <- nodesfilteredsubsists %>%
  group_by(IBD, cluster) %>%
  mutate(ranks = as.numeric(factor(individual)))

nb.cols <- max(nodesfilteredsubsists$ranks[nodesfilteredsubsists$IBD==0.9])
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

gg <- ggplot(edgesfilteredsubsists[edgesfilteredsubsists$IBD==0.9,])+
  geom_rect(data=seasons, aes(xmin=-Inf, xmax=Inf, 
                              ymin=as.POSIXct(date1), ymax=as.POSIXct(date2), 
                              fill=toupper(type)), alpha=0.1)+
  scale_fill_manual("Season",values=c("gold2", "darkblue"))+
  geom_segment(aes(x=groups, xend=groups,
                   y=start, yend=end))+
  geom_point(data = nodesfilteredsubsists[nodesfilteredsubsists$IBD==0.9,],
             aes(x = paste0(IBD,cluster), y = date, color = factor(ranks)),
             show.legend = F, alpha = 0.8, size = 1,
             position = position_jitter(height = 0, width = 0.35, seed = 1234))+
  scale_color_manual(values = mycolors)+
  scale_y_continuous(trans = rev_date, 
                     breaks=unique(as.POSIXct(nodesfiltered$date)), labels = unique(as.Date(nodesfiltered$date)))+
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Sample cluster")+
  ylab("Date")+
  theme(legend.position = NULL,
        axis.text.x = element_text(angle =  90, size=5), 
        panel.grid = element_blank())

gg
ggsave("out/doi-pop-09.png", width = 8, height = 4)

######################################################

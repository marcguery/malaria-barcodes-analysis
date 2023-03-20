######################READING######################
mccoil.snp <- read.table("out/SNP-mccoil_output.txt_summary.txt", h = T)[,-c(1)]
mccoil.wgs <- read.table("out/WGS-mccoil_output.txt_summary.txt", h = T)[,-c(1)]
fws <- read.csv("out/fws-19262-5-4000-persample.csv")
fws.af <- read.table("out/maf.tsv", h = F)
msp <- read.csv("rawdata/MSP-Collins2022.csv", h = T)
msp <- unique(msp[,c("sampleID", "clonenumberbydate")])
msp <- msp[!is.na(msp$sampleID) & msp$clonenumberbydate>0,]
pftest <- read.csv("../../pf-test.csv")

#One sample has sd too high
if (max(mccoil.snp$sd) >1){
  mccoil.snp <- mccoil.snp[mccoil.snp$sd<max(mccoil.snp$sd),]
}
############################################

######################MERGING######################
coi <- merge(mccoil.snp[mccoil.snp$CorP=="C",], 
                mccoil.wgs[mccoil.wgs$CorP=="C",], 
                by = "name", suffixes= c(".SNP", ".WGS"), all = T)

coi <- merge(coi, msp, by.x = "name", by.y = "sampleID", all = TRUE)

coi <- merge(coi, fws, by.x = "name", by.y = "Sample", all = TRUE)


colnames(fws.af) <- c("name", "AF")
af <- mccoil.snp[mccoil.snp$CorP=="P",]
af2 <- mccoil.wgs[mccoil.wgs$CorP=="P",]

af <- data.frame(t(apply(af, c(1),
                         FUN = function(x){
                           if (as.numeric(x[6])+as.numeric(x[7])>=1){
                             x[3] = min(0.5, 1 - as.numeric(x[3]))
                             x[4] = min(0.5, 1 - as.numeric(x[4]))
                             x[6] = min(0.5, 1 -  as.numeric(x[6]))
                             x[7] = min(0.5, 1 -  as.numeric(x[7]))
                           }
                           return(x)
                         })))

af2 <- data.frame(t(apply(af2, c(1),
                          FUN = function(x){
                            if (as.numeric(x[6])+as.numeric(x[7])>=1){
                              x[3] = min(0.5, 1 - as.numeric(x[3]))
                              x[4] = min(0.5, 1 - as.numeric(x[4]))
                              x[6] = min(0.5, 1 -  as.numeric(x[6]))
                              x[7] = min(0.5, 1 -  as.numeric(x[7]))
                            }
                            return(x)
                          })))
write.csv(af[,c(2,3,4)], "out/SNP-MAF-prediction.csv", quote = F, row.names = F)
write.csv(af2[,c(2,3,4)], "out/WGS-MAF-prediction.csv", quote = F, row.names = F)
af <- merge(af,
            af2,
            suffixes= c(".SNP", ".WGS"),
            by = "name", all = T)

af <- af[,-c(which(colnames(af)%in%c("CorP.SNP", "CorP.WGS")))]

af <- merge(af, fws.af, by = "name", all = T)
summary(af$AF - as.numeric(af$median.SNP))
summary(af$AF - as.numeric(af$median.WGS))
############################################

gg <- ggplot(coi)+
  geom_tile(data = data.frame(x = unlist(lapply(1:round(max(coi[,"mean.SNP"], na.rm = T)),
                                                FUN = function(x){
                                                  x+c(0,0,0)})), 
                              y = rep(c(1:round(max(coi[,"mean.WGS"], na.rm = T))),2)),
            aes(x = factor(x), y = factor(y),
                fill = x==y), 
            na.rm = T, color = NA, show.legend = F,
            alpha = 0.5, linetype = 2, size = 1)+
  geom_vline(xintercept = 1.5, linetype = 2, size = 1)+
  geom_hline(yintercept = 1.5, linetype = 2, size = 1)+
  geom_hline(yintercept = 2.5, linetype = 2, size = 1)+
  geom_point(aes(x = round(mean.SNP), y = round(mean.WGS)), na.rm = T, 
             position = position_auto(seed = 456), size = 2, alpha =0.8)+
  xlab("Barcode estimated COI (McCOIL)")+
  ylab("Genome estimated COI (McCOIL)")+
  scale_fill_manual(breaks = c("TRUE", "FALSE"),
                    values = c("green", "red"))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(panel.grid = element_blank(),
        text = element_text(size = 18))
gg
ggsave("out/mccoilSNP-mccoilWGS.png", width = 5, height = 5, units = "in", dpi = 150)

gg <- ggplot(coi)+
  geom_tile(data = data.frame(x = unlist(lapply(1:max(coi[,"clonenumberbydate"], na.rm = T),
                                                      FUN = function(x){
                                                        x+c(rep(0,2))})), 
                              y = rep(c(1:round(max(coi[,"mean.SNP"], na.rm = T))),2)),
            aes(x = factor(x), y = factor(y),
                fill = x==y), 
            na.rm = T, color = NA, show.legend = F,
            alpha = 0.5, linetype = 2, size = 1)+
  geom_vline(xintercept = -0.5+c(2:round(max(coi[,"clonenumberbydate"], na.rm = T))), 
             linetype = 2, size = 1)+
  geom_hline(yintercept = -0.5+c(2:round(max(coi[,"mean.SNP"], na.rm = T))), 
             linetype = 2, size = 1)+
  geom_point(aes(x = round(clonenumberbydate), y = round(mean.SNP)), na.rm = T, 
             position = position_auto(seed = 123), size = 2, alpha =0.8)+
  xlab("MSP2 estimated COI")+
  ylab("Barcode estimated COI (McCOIL)")+
  scale_fill_manual(breaks = c("TRUE", "FALSE"),
                    values = c("green", "red"))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(panel.grid = element_blank(),
        text = element_text(size = 18))
gg
ggsave("out/mccoilSNP-MSP2.png", width = 15, height = 5, units = "in", dpi = 150)

gg <- ggplot(coi)+
  geom_tile(data = data.frame(x = unlist(lapply(1:max(coi[,"clonenumberbydate"], na.rm = T),
                                                FUN = function(x){
                                                  x+c(0,0,0)})), 
                              y = rep(c(1:round(max(coi[,"mean.WGS"], na.rm = T))),2)),
            aes(x = factor(x), y = factor(y),
                fill = x==y), 
            na.rm = T, color = NA, show.legend = F,
            alpha = 0.5, linetype = 2, size = 1)+
  geom_vline(xintercept = -0.5+c(2:round(max(coi[,"clonenumberbydate"], na.rm = T))), 
             linetype = 2, size = 1)+
  geom_hline(yintercept = -0.5+c(2:round(max(coi[!is.na(coi$Fws),"mean.WGS"], na.rm = T))), 
             linetype = 2, size = 1)+
  geom_point(aes(x = round(clonenumberbydate), y = round(mean.WGS)), na.rm = T, 
             position = position_auto(seed = 123), size = 2, alpha =0.8)+
  xlab("MSP2 estimated COI")+
  ylab("Genome estimated COI (McCOIL)")+
  scale_fill_manual(breaks = c("TRUE", "FALSE"),
                    values = c("green", "red"))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(panel.grid = element_blank(),
        text = element_text(size = 18))
gg
ggsave("out/mccoilWGS-MSP2.png", width = 15, height = 5, units = "in", dpi = 150)

coi$clonenumberbydate[coi$clonenumberbydate>2] <- 3
############################################

######################METHOD COMPARISON######################

gg <- ggplot(coi[!is.na(coi$mean.SNP) & !is.na(coi$Fws),])+
  geom_blank(aes(x = factor(mean.SNP),
                 y = Fws))+
  geom_vline(xintercept = 1.5, linetype = 2, size = 1)+
  geom_hline(yintercept = 0.95, color = "blue", linetype = 3)+
  geom_point(aes(x = factor(mean.SNP),
                 y = Fws,
                  color = (Fws > 0.95 & mean.SNP == 1 | Fws < 0.95 & mean.SNP > 1)), 
             position = position_auto(seed = 123), alpha = 0.8, show.legend = F)+
  scale_color_manual(breaks = c("TRUE", "FALSE"),
                    values = c("green", "red"))+
  scale_y_continuous(n.breaks = 10, limits = c(0.3,1.01),
                     expand = c(0,0))+
  xlab("Barcode estimated COI (McCOIL)")+
  ylab("Fws")+
  theme(panel.grid.major.x = element_blank(),
        text = element_text(size = 18))
gg
ggsave("out/mccoilSNP-Fws.png", width = 5, height = 5, units = "in", dpi = 150)

gg <- ggplot(coi[!is.na(coi$mean.WGS) & !is.na(coi$Fws),])+
  geom_blank(aes(x = factor(mean.WGS),
                 y = Fws))+
  geom_vline(xintercept = 1.5, linetype = 2, size = 1)+
  geom_vline(xintercept = 2.5, linetype = 2, size = 1)+
  geom_hline(yintercept = 0.95, color = "blue", linetype = 3)+
  geom_point(aes(x = factor(mean.WGS),
                 y = Fws,
                 color = (Fws > 0.95 & mean.WGS == 1 | Fws < 0.95 & mean.WGS > 1)), 
             position = position_auto(seed = 123), alpha = 0.8, show.legend = F)+
  scale_color_manual(breaks = c("TRUE", "FALSE"),
                     values = c("green", "red"))+
  scale_y_continuous(n.breaks = 10, limits = c(0.3,1.01),
                     expand = c(0,0))+
  xlab("Genome estimated COI (McCOIL)")+
  ylab("Fws")+
  theme(panel.grid.major.x = element_blank(),
        text = element_text(size = 18))
gg
ggsave("out/mccoilWGS-Fws.png", width = 8, height = 5, units = "in", dpi = 150)

gg <- ggplot(coi[!is.na(coi$clonenumberbydate) & !is.na(coi$Fws),])+
  geom_blank(aes(x = factor(clonenumberbydate),
                 y = Fws))+
  geom_vline(xintercept = c(1.5,2.5), linetype = 2, size = 1)+
  geom_hline(yintercept = 0.95, color = "blue", linetype = 3)+
  geom_point(aes(x = factor(clonenumberbydate),
                 y = Fws,
                 group = factor(clonenumberbydate),
                 color = (Fws > 0.95 & clonenumberbydate == 1 | Fws < 0.95 & clonenumberbydate > 1)), 
             position = position_auto(seed = 123), 
             alpha = 0.8, show.legend = F)+
  scale_color_manual(breaks = c("TRUE", "FALSE"),
                     values = c("green", "red"))+
  scale_x_discrete(breaks = c(1,2,3),
                   labels = c(1,2,"3+"))+
  scale_y_continuous(n.breaks = 10, limits = c(0.3,1.01),
                     expand = c(0,0))+
  xlab("MSP2 estimated COI")+
  ylab("Fws")+
  theme(panel.grid.major.x = element_blank(),
        text = element_text(size = 18))
gg
ggsave("out/MSP2-Fws.png", width = 8, height = 5, units = "in", dpi = 150)

coi$date <- paste0(sub("\\S+_", "", coi$name),"01")

coiprop <- coi%>%
  group_by(date)%>%
  summarise(mccoilbarcodecoi = round(length(which(mean.SNP>1))/length(which(!is.na(mean.SNP))),2),
            mccoilbarcodenum = length(which(!is.na(mean.SNP))),
            mccoilwgscoi = round(length(which(mean.WGS>1))/length(which(!is.na(mean.WGS))),2),
            mccoilwgsnum = length(which(!is.na(mean.WGS))),
            fwscoi = round(length(which(Fws<0.95))/length(which(!is.na(Fws))),2),
            fwsnum = length(which(!is.na(Fws))),
            mspcoi = round(length(which(clonenumberbydate>1))/length(which(!is.na(clonenumberbydate))),2),
            mspnum = length(which(!is.na(clonenumberbydate))),
  )

coiprop.melted <- melt(coiprop, id.vars = "date", measure.vars = seq(2,ncol(coiprop),2))
coiprop.melted2 <- melt(coiprop, id.vars = "date", measure.vars = seq(3,ncol(coiprop),2))
coiprop.melted$variable <- sub("coi$", "", coiprop.melted$variable)
coiprop.melted2$variable <- sub("num$", "", coiprop.melted2$variable)
colnames(coiprop.melted) <- c("Date", "Method", "COI")
colnames(coiprop.melted2) <- c("Date", "Method", "Numobs")
coiprop.melted <- merge(coiprop.melted, coiprop.melted2)
coiprop.melted <- coiprop.melted[with(coiprop.melted, order(Method, Date)),]
rm(coiprop.melted2)

coiprop.melted$Date <- as.Date(coiprop.melted$Date, format = "%y%m%d")

#dry (january to july) wet (beginning august to end december)
wetseason <- c(8,12)
wetseason <- as.character(wetseason)

wetseason[nchar(wetseason)==1] <- paste0("0", wetseason[nchar(wetseason)==1])

years <- sort(as.numeric(unique(format(coiprop.melted$Date,"%Y"))))
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
mindate <- min(coiprop.melted$Date)
maxdate <- max(coiprop.melted$Date)
#Time zone warning but it is ok
seasons <- seasons[seasons$date2>mindate & seasons$date1 < maxdate,]
seasons$date2[which.max(seasons$date2)] <- maxdate+15
seasons$date1[which.min(seasons$date1)] <- mindate-15

axisshift <- 0.08
gg <- ggplot(data=coiprop.melted[coiprop.melted$Numobs>=5 & coiprop.melted$Method!="msp",])+
  geom_rect(data=seasons, aes(ymin=-Inf, ymax=Inf, 
                              xmin=date1, xmax=date2, 
                              fill=toupper(type)), alpha=0.5)+
  scale_fill_manual("Season",values=c("gold2", "darkblue"))+
  geom_point(aes(x = Date, y = COI, color = Method))+
  geom_path(aes(x = Date, y = COI, color = Method, group = Method))+
  scale_color_brewer("COI estimation", type = "qual", palette = "Set1",
                     breaks = c("fws", "mccoilbarcode", "mccoilwgs"),
                     labels = c("Fws", "Barcode (McCOIL)", "Genome (McCOIL)"))+
  scale_x_date("", expand = c(0,0), date_breaks = "2 months", date_minor_breaks = "1 month",
               date_labels = "%b")+
  scale_y_continuous(n.breaks = 10)+
  ylab("Proportion of mixed infections")+
  annotate("text",
           label = "2014",
           x=as.Date("2014-12-7"),
           y=-axisshift, 
           size = 3) +
  annotate("segment",
           x=as.Date("2014-11-18"),
           xend=as.Date("2014-12-30"),
           y=-axisshift+0.02,
           yend=-axisshift+0.02,
           color ="grey60") +
  annotate("text",
           label = "2015",
           x=as.Date("2015-06-22"),
           y=-axisshift, 
           size = 3) +
  annotate("segment",
           x=as.Date("2015-01-03"),
           xend=as.Date("2015-12-30"),
           y=-axisshift+0.02,
           yend=-axisshift+0.02,
           color ="grey60") +
  annotate("text",
           label = "2016",
           x=as.Date("2016-06-22"),
           y=-axisshift, 
           size = 3) +
  annotate("segment",
           x=as.Date("2016-01-03"),
           xend=as.Date("2016-12-30"),
           y=-axisshift+0.02,
           yend=-axisshift+0.02,
           color ="grey60") +
  annotate("text",
           label = "2017",
           x=as.Date("2017-03-15"),
           y=-axisshift, 
           size = 3) +
  annotate("segment",
           x=as.Date("2017-01-03"),
           xend=as.Date("2017-05-15"),
           y=-axisshift+0.02,
           yend=-axisshift+0.02,
           color ="grey60") +
  coord_cartesian(ylim=c(-0.01,0.66), clip="off", expand = 0)+
  theme(plot.margin = unit(c(1,1,2,1), "lines"),
        panel.grid = element_line(color = "white"),
        panel.grid.major.x = element_blank())
gg
ggsave("out/coiprop-methods.png", width = 12, height = 4, units = "in", dpi = 150)
#Estimates of proportion of monoclonal infections from:
#mcCOIL on consensus barcode 101 SNPs
round(length(which(coi$mean.SNP==1))/nrow(coi[!is.na(coi$mean.SNP),]),2)
#mcCOIL on WGS 1k SNPs
round(length(which(coi$mean.WGS==1))/nrow(coi[!is.na(coi$mean.WGS),]),2)
#FWS on WGS 19k SNPs
round(length(which(coi$Fws>=0.95))/nrow(coi[!is.na(coi$Fws),]),2)
#MSP2 band lengths
round(length(which(coi$clonenumberbydate>1))/nrow(coi[!is.na(coi$clonenumberbydate),]),2)

#SNP COI vs WGS COI
length(which(coi$mean.SNP==coi$mean.WGS))/length(which(!is.na(coi$mean.SNP==coi$mean.WGS)))
length(which(coi$mean.SNP==1 & coi$mean.WGS>1))/length(which(coi$mean.WGS>1))

#SNP COI and WGS COI vs Fws
length(which(coi$mean.SNP==1 & coi$Fws < 0.95))/length(which(coi$mean.SNP==1))
length(which(coi$mean.WGS==1 & coi$Fws < 0.95))/length(which(coi$mean.WGS==1))

length(which(coi$mean.SNP==1 & coi$Fws < 0.9))/length(which(coi$mean.SNP==1))
length(which(coi$mean.WGS==1 & coi$Fws < 0.9))/length(which(coi$mean.WGS==1))
############################################

####################AF COMPARISON########################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(3)

gg <- ggplot(af)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_dotplot(aes(x = 0,
                   y = -AF + as.numeric(median.SNP)),
               binaxis = "y", stackratio = 1.2, stackdir = "center",
               dotsize = 0.5, binwidth = 0.01,fill = cols[3])+
  geom_dotplot(aes(x = 1,
                   y = -AF + as.numeric(median.WGS)),
               binaxis = "y", stackratio = 1.2, stackdir = "center",
               dotsize = 0.5, binwidth = 0.01,fill = cols[2])+
  scale_y_continuous(n.breaks = 10)+
  scale_x_continuous(breaks = c(0, 1),
                     labels = c("Barcode", "Genome"))+
  xlab("McCOIL COI estimation")+
  ylab("Within-sample minor allele frequency prediction error")+
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 12))

gg
ggsave("out/mccoil-AF-error.png", width = 8, height = 4, units = "in", dpi = 150)

gg <- ggplot(af)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_boxplot(aes(x = -0.0075 + round(AF*20)/20,
                   y = -AF + as.numeric(median.SNP),
                   group = -0.0075 + round(AF*20)/20), fill = cols[3], width = 0.01)+
  geom_boxplot(aes(x = 0.0075 + round(AF*20)/20,
                   y = -AF + as.numeric(median.WGS),
                   group = 0.0075 + round(AF*20)/20), fill = cols[2], width = 0.01)+
  scale_x_continuous(n.breaks = 10)+
  scale_y_continuous(n.breaks = 10, limits = c(-0.2,0.2), expand = c(0,0))+
  ylab("Within-sample minor allele frequency prediction error")+
  xlab("True within-sample minor allele frequency")+
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 12))
gg
ggsave("out/mccoil-AF-error-binned.png", width = 12, height = 4, units = "in", dpi = 150)

gg <- ggplot(af)+
  geom_segment(data= data.frame(1),
               x = 0, y = 0,
               xend = 1, yend = 1, color = "red")+
  geom_point(aes(y = as.numeric(median.WGS), x = as.numeric(AF)))+
  scale_x_continuous(n.breaks = 10, limits = c(-0.01,0.51), expand = c(0,0))+
  scale_y_continuous(n.breaks = 10, limits = c(-0.01,0.51), expand = c(0,0))+
  ylab("Within-sample minor allele frequency prediction")+
  xlab("True within-sample minor allele frequency")+
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 12))
gg
ggsave("out/mccoil-AF-comp.png", width = 12, height = 4, units = "in", dpi = 150)

############################################

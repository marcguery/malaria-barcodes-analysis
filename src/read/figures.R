#######################SAMPLING STATISTICS######################
#Malaria cases
pftest <- read.csv("../../pf-test.csv")
pftest.prev <- pftest%>%
  group_by(date)%>%
  summarise(prevalence = length(which(infectivity=="POS"))/length(infectivity))

malcases <- read.table("../../Gambia-cases.tsv", h = T, sep = "\t")

gg <- ggplot(malcases)+
  geom_line(aes(x = as.Date(paste0(Year, "-01-01")),
                y = Point.cases/Population),
            color = "blue")+
  geom_line(data = pftest.prev[pftest.prev$date < "2017-01-01",],
            aes(x = as.Date(date),
                y =prevalence),
            color = "red")+
  scale_x_date(expand = c(0.01,0),
               breaks = as.Date("1999-01-01")+years(seq(1,21,1)), 
               labels = c(1999)+seq(1,21,1))+
  scale_y_continuous(n.breaks = 10)+
  xlab("Year")+
  ylab("Case prevalence")
gg
ggsave("out/malaria-prevalence.png", width = 8, height = 4, dpi = 100)

#Blood samples
gg <- ggplot(pftest)+
  geom_rect(data=seasons, aes(xmin=date1, xmax=date2, 
                              ymin=0, ymax=Inf, 
                              fill=type), alpha=0.6)+
  scale_fill_manual("", labels=c("Dry season", "Wet season"), 
                    values=c("gold2", "darkblue"))+
  new_scale("fill")+
  geom_vline(data=yearspan,
             aes(xintercept = year), size=1, linetype=2)+
  geom_bar(aes(x = as.Date(date), 
               group = factor(paste0(infectivity, Cohort),
                              levels = rev(c("POSFALSE", "POSTRUE", "NEGFALSE", "NEGTRUE")),
                              labels = rev(c("pos", "poscoh",  "neg",  "neg"))),
               fill = factor(paste0(infectivity, Cohort),
                             levels = rev(c("POSFALSE", "POSTRUE", "NEGFALSE", "NEGTRUE")),
                             labels = rev(c("pos", "poscoh",  "neg",  "neg")))), 
           color = "black")+
  geom_text(data=yearspan,
            aes(x = year+months(manualposition), y = 980, label=year(year)), 
            size=12, fontface="bold")+
  scale_x_date(expand = c(0,0), 
               breaks = as.Date("2014-09-01")+months(seq(2,33,2)), 
               labels = c("Nov", 
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May"), )+
  ylab("Blood samples")+
  xlab("Date")+
  scale_fill_manual("",
                    breaks = c("pos", "poscoh", "neg"),
                    labels=c("Pf+", "Pf+ (Cohort)",  "Neg"), 
                    values=c("red", "orange", "grey60"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1050), n.breaks = 10)+
  theme(text=element_text(size=30),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top")
gg
ggsave("out/study-bloodsamples.png", width = 22, height = 7)
#Sequencing
wgsdata <- read.csv("../../19262-WGS-SNPS.csv")[1]
barcodedata <- read.csv("../../101-genotyped-SNPS.csv")[1]

wgsdata$date <- as.Date(paste(sub("\\S+_", "", wgsdata$ID),"01"), format = "%y%m%d")
wgsdata$type <- "genome"
barcodedata$type <- "barcode"
barcodedata$date <- as.Date(paste(sub("\\S+_", "", barcodedata$ID),"01"), format = "%y%m%d")
malariagendata <- merge(barcodedata, wgsdata, by = "ID", all = T)
malariagendata$date <- pmin(malariagendata$date.x, malariagendata$date.y, na.rm = T)
malariagendata$type <- pmin(malariagendata$type.x, malariagendata$type.y)
malariagendata$type[!is.na(malariagendata$type)] <- "both"
malariagendata$type[is.na(malariagendata$type)] <- pmin(malariagendata$type.x[is.na(malariagendata$type)], 
                                                        malariagendata$type.y[is.na(malariagendata$type)], na.rm = T)
malariagendata <- malariagendata[,c("ID", "date", "type")]
gg <- ggplot(malariagendata)+
  geom_rect(data=seasons, aes(xmin=date1, xmax=date2, 
                              ymin=0, ymax=Inf, 
                              fill=type), alpha=0.6)+
  scale_fill_manual("", labels=c("Dry season", "Wet season"), 
                    values=c("gold2", "darkblue"))+
  new_scale("fill")+
  geom_vline(data=yearspan,
             aes(xintercept = year), size=1, linetype=2)+
  geom_bar(aes(x = as.Date(date), 
               group = factor(type, levels = c("barcode", "genome", "both")),
               fill = type), color = "black")+
  geom_text(data=yearspan,
            aes(x = year+months(manualposition), y = 70, label=year(year)), 
            size=12, fontface="bold")+
  scale_x_date(expand = c(0,0), 
               breaks = as.Date("2014-09-01")+months(seq(2,33,2)), 
               labels = c("Nov", 
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May"), )+
  ylab("Sequences")+
  xlab("Date")+
  scale_fill_manual("",
                    breaks = c("barcode", "genome", "both"),
                    labels=c("Barcode", "Genome", "Both"), 
                    values=c("green", "blue", "red"))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,75), n.breaks = 10,
                     labels = function(x){sprintf("%06s", x)})+
  theme(text=element_text(size=30),
        legend.position = "top",
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank())
gg
ggsave("out/study-malariagen.png", width = 22, height = 7)

gg <- ggplot(meta.goodsnps)+
  
  geom_rect(data=seasons, aes(xmin=date1, xmax=date2, 
                              ymin=0, ymax=Inf, 
                              fill=type), alpha=0.6)+
  geom_vline(data=yearspan,
             aes(xintercept = year), size=1, linetype=2)+
  geom_text(data=yearspan,
            aes(x = year+months(manualposition), y = 70, label=year(year)), 
            size=12, fontface="bold")+
  geom_label(data=data.frame(1), x = as.Date("2017-05-01"), y=30, 
             label="Barcodes", size=8, color="green4", fontface="bold")+
  geom_label(data=data.frame(1), x = as.Date("2017-05-10"), y = 5, 
             label="Genomes", size=8, color="red3", fontface="bold")+
  stat_count(aes(x=Date, color=Study), 
             geom="point", position = "identity", size=2, show.legend = F)+
  stat_count(aes(x=Date, color=Study), 
             geom="line", position = "identity", size=1.5, show.legend = F)+
  scale_x_date(expand = c(0,0), 
               breaks = as.Date("2014-09-01")+months(seq(2,33,2)), 
               labels = c("Nov", 
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May"), )+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,75), n.breaks = 10,
                     labels = function(x){sprintf("%06s", x)})+
  scale_color_manual("Data", labels=c("Genomes", "Barcodes"), values=c("green4", "red3"))+
  scale_fill_manual("", labels=c("Dry season", "Wet season"), 
                    values=c("gold2", "darkblue"))+
  ylab("Samples")+
  xlab("Date")+
  theme(text=element_text(size=30),
        legend.position = "top",
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank())
gg
ggsave("out/study-summary.png", width = 22, height = 7)

meta.dedupl.nocoh <- meta.dedupl[meta.dedupl$Date<"2017-01-01",]
gg <- ggplot()+
  
  geom_rect(data=seasons, aes(xmin=date1, xmax=date2, 
                              ymin=0.5, ymax=Inf, 
                              fill=type), alpha=0.6)+
  geom_vline(data=yearspan,
             aes(xintercept = year), size=1, linetype=2)+
  geom_text(data=yearspan,
            aes(x = year+months(manualposition), y = 70, label=year(year)), 
            size=12, fontface="bold")+
  geom_point(data = data.frame(1:length(unique(meta.dedupl.nocoh$Date))),
             x=unique(meta.dedupl.nocoh$Date), 
             y = 37.5, size=5, pch=21, fill="black")+
  geom_curve(data = data.frame(1),
             x=as.Date("2014-11-01"), xend=as.Date("2015-05-01"), 
             y=45, yend=45, curvature = -0.5, 
             size=3, color="darkgreen",
             arrow = arrow(length=unit(0.8,"cm"), ends="both", type = "closed"))+
  geom_curve(data = data.frame(1),
             x=as.Date("2015-05-01"), xend=as.Date("2015-11-01"), 
             y=30, yend=30, curvature = 0.5, 
             size=3, color="darkblue",
             arrow = arrow(length=unit(0.8,"cm"), ends="both", type = "closed"))+
  geom_curve(data = data.frame(1),
             x=as.Date("2015-11-01"), xend=as.Date("2016-05-01"), 
             y=45, yend=45, curvature = -0.5, 
             size=3, color="darkgreen",
             arrow = arrow(length=unit(0.8,"cm"), ends="both", type = "closed"))+
  geom_curve(data = data.frame(1),
             x=as.Date("2016-05-01"), xend=as.Date("2016-11-01"), 
             y=30, yend=30, curvature = 0.5, 
             size=3, color="darkblue",
             arrow = arrow(length=unit(0.8,"cm"), ends="both", type = "closed"))+
  scale_x_date(expand = c(0,0), 
               breaks = as.Date("2014-09-01")+months(seq(2,33,2)), 
               labels = c("Nov", 
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May"), )+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,75), n.breaks = 10,
                     labels = function(x){sprintf("%06s", x)})+
  scale_fill_manual("", labels=c("Dry season", "Wet season"), values=c("gold2", "darkblue"))+
  ylab("")+
  xlab("Date")+
  theme(text=element_text(size=30),
        axis.text.y = element_text(color = NA),
        legend.position= "top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.y = element_blank())
gg
ggsave("out/wetdry-summary.png", width = 22, height = 7)
############################################
#######################SNP LOCATION######################

gg <- ggplot(data = snps.dedupl)+
  geom_rect(data=chs,
            aes(xmin = 0, xmax=length, 
                ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
            color="grey70", fill=NA)+
  geom_text(data=chs,
            aes(x = length+100000, y=as.numeric(chr)*2.75, 
                label=paste(Usedsnp, "/", Tsnp)), 
            color="grey50", size=3)+
  geom_text(data=chs,
            aes(x = 0-200000, y=as.numeric(chr)*2.75, 
                label=chr), 
            color="grey50", size=3)+
  geom_segment(data = snps.dedupl[seq(1, nrow(snps.dedupl), 2),],
               aes(x = pos, xend=pos, 
                   y=as.numeric(chr)*2.75-0.5, yend=as.numeric(chr)*2.75-0.1, 
                   color=meanComp, alpha = avail), 
               size=0.4)+
  geom_segment(data = snps.dedupl[seq(2, nrow(snps.dedupl), 2),],
               aes(x = pos, xend=pos, 
                   y=as.numeric(chr)*2.75+0.5, yend=as.numeric(chr)*2.75+0.1, 
                   color=meanComp, alpha = avail), 
               size=0.4)+
  scale_color_steps2(midpoint = 0.8, high="springgreen4", mid = "orange2", low="red")+
  theme(plot.background = element_blank(),panel.grid = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
gg

ggsave("out/snp-concordance-dedupl.png", width =15, height = 8)

gg <- ggplot(data = snps.goodsnps)+
  geom_rect(data=chs,
            aes(xmin = 0, xmax=length, 
                ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
            color="grey70", fill=NA)+
  geom_text(data=chs,
            aes(x = length+100000, y=as.numeric(chr)*2.75, 
                label=paste(Usedsnp, "/", Tsnp)), 
            color="grey50", size=3)+
  geom_text(data=chs,
            aes(x = 0-200000, y=as.numeric(chr)*2.75, 
                label=chr), 
            color="grey50", size=3)+
  geom_segment(data = snps.goodsnps[seq(1, nrow(snps.goodsnps), 2),],
               aes(x = pos, xend=pos, 
                   y=as.numeric(chr)*2.75-0.5, yend=as.numeric(chr)*2.75-0.1, 
                   color=meanComp, alpha = avail), 
               size=0.4)+
  geom_segment(data = snps.goodsnps[seq(2, nrow(snps.goodsnps), 2),],
               aes(x = pos, xend=pos, 
                   y=as.numeric(chr)*2.75+0.5, yend=as.numeric(chr)*2.75+0.1, 
                   color=meanComp, alpha = avail), 
               size=0.4)+
  scale_color_steps2(midpoint = 0.8, high="springgreen4", mid = "orange2", low="red")+
  theme(plot.background = element_blank(),panel.grid = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
gg
ggsave("out/snp-concordance-goodsnps.png", width =15, height = 8)

gg <- ggplot(data = snps.cons)+
  geom_rect(data=chs,
            aes(xmin = 0, xmax=length, 
                ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
            color="grey70", fill=NA)+
  geom_text(data=chs,
            aes(x = length+100000, y=as.numeric(chr)*2.75, 
                label=paste(Usedsnp, "/", Tsnp)), 
            color="grey50", size=3)+
  geom_text(data=chs,
            aes(x = 0-200000, y=as.numeric(chr)*2.75, 
                label=chr), 
            color="grey50", size=3)+
  geom_segment(data = snps.cons[seq(1, nrow(snps.cons), 2),],
               aes(x = pos, xend=pos, 
                   y=as.numeric(chr)*2.75-0.5, yend=as.numeric(chr)*2.75-0.1, 
                   color=meanComp, alpha = avail), 
               size=0.4)+
  geom_segment(data = snps.cons[seq(2, nrow(snps.cons), 2),],
               aes(x = pos, xend=pos, 
                   y=as.numeric(chr)*2.75+0.5, yend=as.numeric(chr)*2.75+0.1, 
                   color=meanComp, alpha = avail), 
               size=0.4)+
  scale_color_steps2(midpoint = 0.8, high="springgreen4", mid = "orange2", low="red", na.value = "black")+
  theme(plot.background = element_blank(),panel.grid = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
gg
ggsave("out/snp-concordance-consensus.png", width =15, height = 8)
############################################

######################FILTERING 21 BAD GENOTYPED SNPS######################
#BEFORE REMOVAL

#These 21 bad SNPs correlate well with WGS SNP before May 2016
gg <- ggplot(data=subset(bcodescores.dedupl.melted, date<=1605))+
  geom_tile(aes(x=factor(Var1, levels=orderedsnps), y=factor(Var2, levels=orderedbarcodes),
                fill=factor(value, levels=c("-2","-1","0","1","2"))))+
  coord_fixed(expand = T)+
  scale_fill_manual(labels=c("Match", "Mismatch", "N-N", "N-ATGC", "X-ATGCNX"),
                    values=c("1" = "green4", "-1" = "red3", 
                             "2" = "cyan3", "-2" = "yellow2",
                             "0" = "grey60"
                             ),
                    drop = FALSE,
                    name = "")+
  xlab("SNPs")+
  ylab("Barcodes")+
  theme(axis.text.x = element_text(angle = 90))

gg
ggsave("out/before-removal-before0516.png", plot = gg,
       width = 14, height = 10)

#These 21 bad SNPs correlate badly with WGS SNP after May 2016
gg <- ggplot(data=subset(bcodescores.dedupl.melted, date>1605))+
  geom_tile(aes(x=factor(Var1, levels=orderedsnps), y=factor(Var2, levels=orderedbarcodes),
                fill=factor(value, levels=c("-2","-1","0","1","2"))))+
  coord_fixed(expand = T)+
  scale_fill_manual(labels=c("Match", "Mismatch", "N-N", "N-ATGC", "X-ATGCNX"),
                    values=c("1" = "green4", "-1" = "red3", 
                             "2" = "cyan3", "-2" = "yellow2",
                             "0" = "grey60"
                    ),
                    drop = FALSE,
                    name = "")+
  xlab("SNPs")+
  ylab("Barcodes")+
  theme(axis.text.x = element_text(angle = 90))

gg
ggsave("out/before-removal-after0516.png", plot = gg,
       width = 15, height = 22)

#AFTER REMOVAL
#Before May 2016
gg <- ggplot(data=subset(bcodescores.goodsnps.melted, date<=1605))+
  geom_tile(aes(x=factor(Var1, levels=orderedsnps), y=factor(Var2, levels=orderedbarcodes),
                fill=factor(value, levels=c("-2","-1","0","1","2"))))+
  coord_fixed(expand = T)+
  scale_fill_manual(labels=c("Match", "Mismatch", "N-N", "N-ATGC", "X-ATGCNX"),
                    values=c("1" = "green4", "-1" = "red3", 
                             "2" = "cyan3", "-2" = "yellow2",
                             "0" = "grey60"
                    ),
                    drop = FALSE,
                    name = "")+
  xlab("SNPs")+
  ylab("Corrected barcodes")+
  theme(axis.text.x = element_text(angle = 90, size=3),
        axis.text.y = element_text(size=3))

gg
#After May 2016
gg <- ggplot(data=subset(bcodescores.goodsnps.melted, date>1605))+
  geom_tile(aes(x=factor(Var1, levels=orderedsnps), y=factor(Var2, levels=orderedbarcodes),
                fill=factor(value, levels=c("-2","-1","0","1","2"))))+
  coord_fixed(expand = T)+
  scale_fill_manual(labels=c("Match", "Mismatch", "N-N", "N-ATGC", "X-ATGCNX"),
                    values=c("1" = "green4", "-1" = "red3", 
                             "2" = "cyan3", "-2" = "yellow2",
                             "0" = "grey60"
                    ),
                    drop = FALSE,
                    name = "")+
  xlab("SNPs")+
  ylab("Corrected barcodes")+
  theme(axis.text.x = element_text(angle = 90, size=3),
        axis.text.y = element_text(size=3))

gg
############################################

###################MAKING CONSENSUS BARCODES###################
#Before May 2016
gg <- ggplot(data=subset(bcodescores.cons.melted, date<=1605))+
  geom_tile(aes(x=factor(Var1, levels=orderedsnps), y=factor(Var2, levels=orderedbarcodes),
                fill=factor(value, levels=c("-2","-1","0","1","2"))))+
  coord_fixed(expand = T)+
  scale_fill_manual(labels=c("Match", "Mismatch", "N-N", "N-ATGC", "X-ATGCNX"),
                    values=c("1" = "green4", "-1" = "red3", 
                             "2" = "cyan3", "-2" = "yellow2",
                             "0" = "grey60"
                    ),
                    drop = FALSE,
                    name = "")+
  xlab("SNPs")+
  ylab("Consensus Barcodes")+
  theme(axis.text.x = element_text(angle = 90))

gg
ggsave("out/after-removal-before0516.png", plot = gg,
       width = 14, height = 10)

#After May 2016
gg <- ggplot(data=subset(bcodescores.cons.melted, date>1605))+
  geom_tile(aes(x=factor(Var1, levels=orderedsnps), y=factor(Var2, levels=orderedbarcodes),
                fill=factor(value, levels=c("-2","-1","0","1","2"))))+
  coord_fixed(expand = T)+
  scale_fill_manual(labels=c("Match", "Mismatch", "N-N", "N-ATGC", "X-ATGCNX"),
                    values=c("1" = "green4", "-1" = "red3", 
                             "2" = "cyan3", "-2" = "yellow2",
                             "0" = "grey60"
                    ),
                    drop = FALSE,
                    name = "")+
  xlab("SNPs")+
  ylab("Consensus Barcodes")+
  theme(axis.text.x = element_text(angle = 90))

gg
ggsave("out/after-removal-after0516.png", plot = gg,
       width = 15, height = 22)
######################################

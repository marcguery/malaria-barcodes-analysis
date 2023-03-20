###########STATS##############
if(removecontinf){
  outfolder="out/nocontinf/"
}else
{
  outfolder="out/"
}

#Number of barcodes involved in the IBD network
length(unique(c(databarcodes$sample2, databarcodes$sample1)))
#Number of IBD comparisons with at least 10 informative sites
length(which(databarcodes$N_informative_sites>=10))

yearspan <- data.frame(as.Date("2014-01-01")+months(c(0,12,24,36)))
colnames(yearspan) <- "year"
yearspan$manualposition <- c(9, 6, 6, 3)

gg <- ggplot(relatednessscore.date)+
  geom_rect(data=seasons, aes(xmin=date1, xmax=date2, 
                              ymin=0, ymax=Inf, 
                              fill=type), alpha=0.6)+
  scale_fill_manual("", labels=c("Dry season", "Wet season"), values=c("gold2", "darkblue"))+
  new_scale("fill")+
  geom_vline(data=yearspan,
             aes(xintercept = year), size=1, linetype=2)+
  geom_col(aes(x = as.Date(date), y = numobs,
               fill = factor(1-round(propshared*10)/10)), 
           color = "black", position = "stack")+
  geom_text(data=yearspan,
            aes(x = year+months(manualposition), y = 61.6, label=year(year)), 
            size=12, fontface="bold")+
  scale_x_date(expand = c(0,0), 
               breaks = as.Date("2014-09-01")+months(seq(2,33,2)), 
               labels = c("Nov", 
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May"), )+
  ylab("Number of barcodes")+
  xlab("Date")+
  scale_fill_brewer("Proportion of unique barcodes:", palette = "Reds", type = "seq")+
  scale_y_continuous(expand = c(0,0),limits = c(0,66), n.breaks = 10,
                     labels = function(x){sprintf("%06s", x)})+
  theme(text=element_text(size=30),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top")
gg
ggsave("out/unique-barcodes.png", width = 22, height = 7)
############################

##############FIGURES##############
gg <- ggplot(edgesfiltered.all.nocoh[!edgesfiltered.all.nocoh$sameIndividual,])+
  stat_summary(aes(x = paste0(sameVillage, sameCompound),
                   y = fract_sites_IBD,
                   fill = paste0(sameVillage, sameCompound),
                   group = paste0(sameVillage, sameCompound)),
               fun = function(x){
                 length(which(x>0.5))/length(x)},
               show.legend = F, geom = "col" ,color = "black")+
  ylab("Proportion of related pairs")+
  scale_fill_brewer(name = "",
                    breaks = c("TRUETRUE", "TRUEFALSE", "FALSEFALSE"),
                    labels = list("TRUETRUE"= "Same household",
                                  "TRUEFALSE" = "Different households (same village)",
                                  "FALSEFALSE" = "Different households (different villages)"),
                    palette = 9, type = "seq")+
  scale_x_discrete(labels = c("Different households (different villages)", 
                              "Different households (same village)", 
                              "Same household"))+
  coord_cartesian(ylim = c(0,0.0425), expand = c(0))+
  theme(text = element_text(size=15),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
gg


ggsave(paste0(outfolder, "compound-similarity.png"), plot = gg,
       width = 12, height = 6)

villagediff <- chisq.test(table(edgesfiltered.all.nocoh$sameVillage[!edgesfiltered.all.nocoh$sameIndividual], 
                         edgesfiltered.all.nocoh$fract_sites_IBD[!edgesfiltered.all.nocoh$sameIndividual] > 0))
villagediff

compounddiff <- chisq.test(table(edgesfiltered.all.nocoh$sameCompound[!edgesfiltered.all.nocoh$sameIndividual], 
                         edgesfiltered.all.nocoh$fract_sites_IBD[!edgesfiltered.all.nocoh$sameIndividual] > 0))
compounddiff

edges.grouped.compound$compound1 <- substr(edges.grouped.compound$commonCompound, 1, 4)
edges.grouped.compound$compound2 <- ifelse(nchar(edges.grouped.compound$commonCompound) == 4,
                                           substr(edges.grouped.compound$commonCompound, 1, 4),
                                           substr(edges.grouped.compound$commonCompound, 5, 8))
gg <- ggplot(edges.grouped.compound)+
  geom_blank(data = data.frame(unique(c(edges.grouped.compound$compound1, 
                                        edges.grouped.compound$compound2))),
             mapping = aes(x = unique(c(edges.grouped.compound$compound1, 
                                        edges.grouped.compound$compound2)),
                           y = unique(c(edges.grouped.compound$compound1, 
                                        edges.grouped.compound$compound2))))+
  geom_segment(data = data.frame(1),
               aes(x = 0, y = 0,
                   xend = Inf,
                   yend = Inf),
               linetype = 2, size = 0.3)+
  geom_point(data = edges.grouped.compound[edges.grouped.compound$horizontalIBD != 0,],
             aes(x = compound1,
                 y = compound2,
                 fill = horizontalIBD/horizontalMax,
                 color = commonVillage),
             pch = 21, show.legend = F)+
  geom_point(data = edges.grouped.compound[edges.grouped.compound$horizontalIBD == 0,],
             aes(x = compound1,
                 y = compound2),
             pch = 21, fill = "grey90", color = "grey90")+
  scale_fill_gradient(low = "white", high = "black")+
  scale_y_discrete()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank())+
  guides(color = guide_legend(ncol = 2))
  
gg

ggsave(paste0(outfolder, "compound-comparison.png"), plot = gg,
       width = 8, height = 8)

##############DATES################
colforbars <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))

gg <- ggplot(edges.grouped.date)+
  geom_rect(data= dfgrouplocation,
            aes(xmin=mingrouplocation, xmax = maxgrouplocation,
                ymin = 0, ymax = meanIBD,
                group = group,
              fill = group), 
           position = "dodge2", color="black", show.legend = F)+
  stat_summary(aes(x = grouplocation, group = group,
                   y = horizontalIBD/horizontalMax), 
               fun.data = mean_se, geom = "linerange", size = 1)+
  geom_point(aes(x=round(dayElapsed/30), 
                 y=horizontalIBD/horizontalMax), 
             size=4, pch=21, fill="black", color="black")+
  geom_point(aes(x=round(dayElapsed/30), 
                 y=horizontalIBD/horizontalMax, 
                 fill=group), 
             size=3, pch=21, show.legend = F)+
  scale_fill_gradientn("Time span (months)", colors=colforbars(16))+
  scale_x_continuous(breaks = seq(0,
                                  max(unique(round(edges.grouped.date$dayElapsed/30)))),
                     labels=seq(0,
                                max(unique(round(edges.grouped.date$dayElapsed/30)))), expand=c(0.02,0))+
  scale_y_continuous(breaks=seq(0,0.09, 0.01), expand = c(0.02,0))+
  ylab("Mean IBD")+
  xlab("Time span (months)")+
  theme(text = element_text(size=15), panel.grid.minor = element_blank())

gg
ggsave(paste0(outfolder, "IBD-decrease-time.png"), gg, width = 12, height = 6)

gg <- ggplot(edges.grouped.datecompound)+
  geom_boxplot(aes(x=group, 
                 y=horizontalIBD/horizontalMax, 
                 fill=paste0(sameVillage, sameCompound),
                 group=factor(paste0(group, sameVillage, sameCompound),
                              levels = as.vector(outer(unique(sort(group)), 
                                                         c("TRUETRUE", "TRUEFALSE", "FALSEFALSE"), 
                                             paste0)))), 
             outlier.shape = NA)+
  geom_point(aes(x=group, 
                 y=horizontalIBD/horizontalMax, 
                 fill=paste0(sameVillage, sameCompound),
                 group=factor(paste0(group, sameVillage, sameCompound),
                              levels =as.vector(outer(unique(sort(group)), 
                                                        c("TRUETRUE", "TRUEFALSE", "FALSEFALSE"), 
                                            paste0)))), 
             pch=21, show.legend = F, position=position_jitterdodge(jitter.width = 0.1, seed = 123))+
  # stat_summary(aes(x = group,
  #                  y = horizontalIBD/horizontalMax,
  #                  group = paste0(sameVillage, sameCompound),
  #                  color = paste0(sameVillage, sameCompound)),
  #              fun = median, geom = "line", position = position_dodge(width = 0.75),
  #              show.legend = F, linetype = 2)+
  stat_summary(aes(x = group,
                   y = horizontalIBD/horizontalMax,
                   color = "red"), 
               fun = mean, geom = "path", linetype = 2)+
  scale_color_manual(name = "", 
                     labels = c("Mean of all sampling locations"),
                     values = c("red"))+
  stat_summary(aes(x = group,
                   y = horizontalIBD/horizontalMax), 
               fun = mean, geom = "point", color = "red")+
  scale_x_continuous(breaks = dfgrouplocation.datecompound$group,
                     labels=paste0(dfgrouplocation.datecompound$mingrouplocation, 
                                   "-",
                                   dfgrouplocation.datecompound$maxgrouplocation))+
  scale_y_continuous(breaks=seq(0,0.20, 0.01), expand = c(0.01,0))+
  scale_fill_brewer(name = "",
                    breaks = c("TRUETRUE", "TRUEFALSE", "FALSEFALSE"),
                    labels = list("TRUETRUE"= "Same household", 
                               "TRUEFALSE" = "Different households (same village)", 
                               "FALSEFALSE" = "Different households (different villages)"),
                    palette = 9, type = "seq")+
  ylab("Proportion of related isolates")+
  xlab("Months between sampled pairs")+
  theme(text = element_text(size=15), panel.grid.minor = element_blank(),
        legend.position=c(.8,.8), legend.key.size = unit(2, "lines"), legend.spacing.x = unit(1, "lines"),
        legend.background = element_rect(fill = NA), legend.text = element_text(size=11),
        legend.box.background = element_rect(fill = "grey92", color ="black"),
        legend.key=element_blank(), legend.title = element_blank(),
        legend.spacing.y = unit(0, "lines"), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(r = 0.1, l = 0.1, t = 0.1, b = 0.1, unit = "lines"))

gg
ggsave(paste0(outfolder, "compound-time.png"), gg, width = 12, height = 6)

rowsthatmatter <- edges.grouped.datecompound$group==1 & edges.grouped.datecompound$sameCompound
samecless2 <- edges.grouped.datecompound$horizontalIBD[rowsthatmatter]/edges.grouped.datecompound$horizontalMax[rowsthatmatter]
rowsthatmatter <- edges.grouped.datecompound$group==1 & edges.grouped.datecompound$sameVillage & !edges.grouped.datecompound$sameCompound
diffcsamevless2 <- edges.grouped.datecompound$horizontalIBD[rowsthatmatter]/edges.grouped.datecompound$horizontalMax[rowsthatmatter]
rowsthatmatter <- edges.grouped.datecompound$group==1 & !edges.grouped.datecompound$sameVillage
diffvless2 <- edges.grouped.datecompound$horizontalIBD[rowsthatmatter]/edges.grouped.datecompound$horizontalMax[rowsthatmatter]

t.test(samecless2, diffcsamevless2, alternative = "two.sided", var.equal = FALSE)
t.test(samecless2, diffvless2, alternative = "two.sided", var.equal = FALSE)
t.test(diffcsamevless2, diffvless2, alternative = "two.sided", var.equal = FALSE)

morethanonyear <- edges.grouped.date[edges.grouped.date$dayElapsed/30>11.5,]
quantile(morethanonyear$horizontalIBD/morethanonyear$horizontalMax, c(0.05, 0.5, 0.95))
mean(morethanonyear$horizontalIBD/morethanonyear$horizontalMax)



gg <- ggplot(edges.grouped.date.1season)+
  geom_blank(aes(color=paste0(commonSeason, cycleElapsed)))+
  geom_violin(aes(x=factor(paste0(commonSeason, cycleElapsed),
                           levels = c("wet0", "drywet0", "dry0", "drywet1")),
                  y=horizontalIBD/horizontalMax))+
  geom_point(aes(x=paste0(commonSeason, cycleElapsed),
                 y=horizontalIBD/horizontalMax, 
                 fill=paste0(commonSeason, cycleElapsed)), 
             pch=21, size=3, position=position_jitter(width=0.3, seed=100, 
                                                      height = 0), show.legend = F)+
  geom_segment(data=drywettests,
               aes(x = to, xend = to,
                   y = c(0.065,0.075,0.07),
                   yend = c(0.0625,0.0725,0.0675),
                   color = to), linetype = 1, show.legend = F)+
  geom_segment(data=drywettests,
               aes(x = from, xend = from,
                   y = c(0.065,0.075,0.07),
                   yend = c(0.0625,0.0725,0.0675),
                   color = from), linetype = 1, show.legend = F)+
  geom_segment(data=drywettests,
               aes(x = from, xend = to,
                   y = c(0.065,0.075,0.07),
                   yend = c(0.065,0.075,0.07)), linetype = 2, show.legend = F)+
  geom_text(data=drywettests,
               aes(x = c(4.25,4.25,4.25),
                   y = c(0.065,0.075,0.07),
                   label = paste0("p < ",scientific(ceiling(pvalues*1000)/1000, 1))))+
  scale_color_brewer(palette=6, type="qual")+
  scale_fill_brewer(palette=6, type="qual")+
  scale_x_discrete(breaks = c("wet0", "drywet0", "dry0", "drywet1"),
                   labels = c("Intraseasonal wet",
                              "Wet to dry", 
                              "Intraseasonal dry",
                              "Dry to wet"))+
  scale_y_continuous(breaks = c(seq(0,1,0.01)),
                     limits = c(-0.002,0.0775), expand = c(0,0))+
  xlab("Seasonal combination")+
  ylab("Proportion of similar barcodes")+
  theme(text=element_text(size=15),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
gg
ggsave(paste0(outfolder, "wetdry-comparison.png"), plot = gg,
       width = 10, height = 6)

seasons.combination.nocoh$groups <- factor(paste0(seasons.combination.nocoh$season1, 
                                                  seasons.combination.nocoh$season2),
                                     levels = c("drydry", "wetdry",
                                                "drywet", "wetwet"))
gg <- ggplot(data = edges.grouped.date)+
  geom_segment(data = seasons.combination.nocoh,
               aes(y = as.POSIXct(date1start), 
                   yend = as.POSIXct(date1end),
                   color = season1,
                   x = as.POSIXct(date2start),
                   xend = as.POSIXct(date2start)),
               show.legend = NA, linetype = 1, size = 1.25)+
  geom_segment(data = seasons.combination.nocoh[seasons.combination.nocoh$date2start>min(seasons.combination.nocoh$date2start),],
               aes(x = as.POSIXct(date1start), 
                   xend = as.POSIXct(date1end),
                   color = season1,
                   y = as.POSIXct(date2start),
                   yend = as.POSIXct(date2start)),
               show.legend = F, linetype = 1, size = 1.25)+
  geom_polygon(data = data.frame(x = c(min(seasons.nocoh$date1), max(seasons.nocoh$date2), max(seasons.nocoh$date2)), 
                                 y = c(min(seasons.nocoh$date1), min(seasons.nocoh$date1), max(seasons.nocoh$date2))),
               aes(x = as.POSIXct(x), y = as.POSIXct(y)),
               fill = "white", color = "white")+
  geom_segment(data = seasons.nocoh,
               aes(x = as.POSIXct(date1), 
                   xend = as.POSIXct(date2),
                   color = type,
                   y = as.POSIXct(max(date2)),
                   yend = as.POSIXct(max(date2))), linetype = 1, 
               show.legend = F, size = 1.25)+
  geom_point(data = seasons.combination.nocoh,
             aes(y = as.POSIXct(date1start),
                 x = as.POSIXct(date2start)),
             show.legend = NA, size = 1.25, color = "white", pch = 15)+
  geom_point(data = seasons.nocoh,
             aes(x = as.POSIXct(date1),
                 y = as.POSIXct(max(date2))),
             show.legend = NA, size = 1.25, color = "white", pch = 15)+
  scale_color_manual("Season", values = c("orange", "steelblue"), labels = c("Dry", "Wet"))+
  new_scale("color")+
  geom_point(data = edges.grouped.date.not1season,
             aes(x = as.POSIXct(datemin),
                 y = as.POSIXct(datemax),
                 size = horizontalMax),
             color = "black", show.legend = F)+
  geom_point(data = edges.grouped.date.1season,
             aes(x = as.POSIXct(datemin),
                 y = as.POSIXct(datemax),
                 size = horizontalMax),
             color = "black", show.legend = F)+
  geom_point(aes(x = as.POSIXct(datemin),
                 y = as.POSIXct(datemax),
                 size = 0.1*horizontalMax^(1.25),
                 color = horizontalIBD/horizontalMax),
             show.legend = NA)+
  scale_color_gradient("Proportion of similar barcodes", low = "white", high = "black", breaks = seq(0,1,0.02))+
  scale_size_continuous("Number of comparisons", breaks = seq(0,5000,500))+
  scale_y_continuous(trans = rev_date, 
                     breaks = unique(as.POSIXct(c(edges.grouped.date.not1season$datemin, edges.grouped.date.not1season$datemax))), 
                     labels = month(unique(as.Date(c(edges.grouped.date.not1season$datemin, 
                                                     edges.grouped.date.not1season$datemax))), label = T),
                     expand = c(0.01,0))+
  scale_x_continuous(trans = ori_date, 
                     breaks = unique(as.POSIXct(c(edges.grouped.date.not1season$datemin, edges.grouped.date.not1season$datemax))), 
                     labels = month(unique(as.Date(c(edges.grouped.date.not1season$datemin, 
                                                     edges.grouped.date.not1season$datemax))), label = T),
                     expand = c(0.01,0))+
  annotate("text",
           label = "2014",
           x=as.POSIXct("2014-10-23"),
           y=as.POSIXct("2017-01-30"), 
           size = 3.5) +
  annotate("segment",
           x=as.POSIXct("2014-08-18"),
           xend=as.POSIXct("2014-12-30"),
           y=as.POSIXct("2017-01-15"),
           yend=as.POSIXct("2017-01-15"),
           color ="grey60") +
  annotate("text",
           label = "2015",
           x=as.POSIXct("2015-06-15"),
           y=as.POSIXct("2017-01-30"), 
           size = 3.5) +
  annotate("segment",
           x=as.POSIXct("2015-01-03"),
           xend=as.POSIXct("2015-12-30"),
           y=as.POSIXct("2017-01-15"),
           yend=as.POSIXct("2017-01-15"),
           color ="grey60") +
  annotate("text",
           label = "2016",
           x=as.POSIXct("2016-06-08"),
           y=as.POSIXct("2017-01-30"), 
           size = 3.5) +
  annotate("segment",
           x=as.POSIXct("2016-01-03"),
           xend=as.POSIXct("2016-12-15"),
           y=as.POSIXct("2017-01-15"),
           yend=as.POSIXct("2017-01-15"),
           color ="grey60") +
  annotate("text",
           label = "2014",
           y=as.POSIXct("2014-10-23"),
           x=as.POSIXct("2014-06-30"), 
           size = 3.5, angle = -90, hjust = 0.5) +
  annotate("segment",
           y=as.POSIXct("2014-08-18"),
           yend=as.POSIXct("2014-12-30"),
           x=as.POSIXct("2014-07-15"),
           xend=as.POSIXct("2014-07-15"),
           color ="grey60") +
  annotate("text",
           label = "2015",
           y=as.POSIXct("2015-06-15"),
           x=as.POSIXct("2014-06-30"), 
           size = 3.5, angle = -90, hjust = 0.5) +
  annotate("segment",
           y=as.POSIXct("2015-01-03"),
           yend=as.POSIXct("2015-12-30"),
           x=as.POSIXct("2014-07-15"),
           xend=as.POSIXct("2014-07-15"),
           color ="grey60") +
  annotate("text",
           label = "2016",
           y=as.POSIXct("2016-06-08"),
           x=as.POSIXct("2014-06-30"), 
           size = 3.5, angle = -90, hjust = 0.5) +
  annotate("segment",
           y=as.POSIXct("2016-01-03"),
           yend=as.POSIXct("2016-12-15"),
           x=as.POSIXct("2014-07-15"),
           xend=as.POSIXct("2014-07-15"),
           color ="grey60") +
  theme_classic()+
  coord_cartesian(ylim = c(as.POSIXct("2016-12-15"), 
                           as.POSIXct("2014-08-15")),
                  xlim = c(as.POSIXct("2014-08-15"),
                           as.POSIXct("2016-12-15")),
                  clip="off")+
  theme(plot.margin = unit(c(1,1,2,2), "lines"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.85, 0.75),
        legend.key.width = unit(0.5, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        axis.text.y = element_text(angle = -90, hjust = 0.5),
        text = element_text(size = 11),
        axis.line = element_blank(),
        legend.background = element_rect(color = "grey20", size = 0.5))
gg
ggsave(paste0(outfolder, "seasons-combinations.png"), plot = gg,
       width = 9, height = 9)
##################################

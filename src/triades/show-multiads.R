
#############Isolate from same individual#############
family.all <- c("P0080842_1612", 
                "P0080842_1703", 
                "P0080842_1704",
                "P0080842_1705")
for (mem in family.all){
  gg <- familymapv2(chrombarcodes, mem, 
              family.all, 
              chromlength)
  ggsave(paste0("out/", mem, "_with_family.png"), 
         plot = gg, width = 8, height = 5)
}
edges.sameind.all <- edges[edges$sameIndividual & (as.Date(edges$date1) >= as.Date("2016-12-01") & as.Date(edges$date2) >= as.Date("2016-12-01")),]
ind.same <- unique(edges.sameind.all$commonIndividual)
for (ind in ind.same){
  edges.sameind <- edges.sameind.all[edges.sameind.all$commonIndividual == ind,]
  edges.sameind$inverted <- as.Date(edges.sameind$date1) > as.Date(edges.sameind$date2)
  edges.sameind.ibd <- edges.sameind[edges.sameind$inverted,c("date1", "date2", "fract_sites_IBD")]
  colnames(edges.sameind.ibd) <- c("date2", "date1", "fract_sites_IBD")
  edges.sameind.ibd <- rbind(edges.sameind.ibd,edges.sameind[!edges.sameind$inverted,c("date1", "date2", "fract_sites_IBD")])
  gg <- ggplot(edges.sameind.ibd)+
    geom_blank(data = data.frame(c(1,1)),
               aes(fill = c(0,1)))+
    geom_point(aes(y = as.Date(date1),
                  x = as.Date(date2),
                  fill = fract_sites_IBD), pch = 21, 
               size = 25)+
    geom_text(aes(y = as.Date(date1),
                  x = as.Date(date2),
                  label = round(fract_sites_IBD*100)/100))+
    scale_fill_gradient2("IBD",low = "red", mid = "green",
                         high = "green", 
                         midpoint = 0.9)+
    scale_x_date(expand = c(0,0),
                 breaks = as.Date("2016-12-01")+months(seq(1,5,1)),
                 labels = function(x){month(x, label = TRUE,abbr = TRUE)},
                 limits = c(as.Date("2016-12-15"),
                            as.Date("2017-05-15")))+
    scale_y_date(expand = c(0,0),
                 breaks = as.Date("2016-11-01")+months(seq(1,5,1)),
                 labels = function(x){month(x, label = TRUE,abbr = TRUE)},
                 limits = c(as.Date("2016-11-15"),
                            as.Date("2017-04-15")))+
    xlab("")+
    ylab("")+
    ggtitle(ind)+
    theme(legend.position = "top",
          panel.grid = element_blank())
  ggsave(paste0("out/", ind, "IBDmatrix.png"), gg, width = 8, height = 8, dpi = 100)
  
}


##########################

#############Same parasitic parents#############

#IBD of parents

gg <- parentsmapv3(chrombarcodes, 
                   "K0141403_1504", 
                   "K0374502_1412", "K0374502_1412",
                   chlen = chromlength, mode = 1)
gg
ggsave(paste0("out/offsprings/K0141403_1504_with_K0374502_1412.png"), 
       plot = gg, width = 8, height = 5)

#

siblings.all <- c("J0060701_1511",
                  "J0262731_1412", 
                  "K0030301_1504",
                  "K0030302_1512",
                  "K0161621_1506",
                  "K0374408_1412")
for (sib in siblings.all){
  gg <- familymapv2(chrombarcodes, sib, siblings.all, 
              chromlength)
  ggsave(paste0("out/siblings/", sib, "_with_siblings.png"), 
         plot = gg, width = 8, height = 5)
}
gg <- familymapv2(chrombarcodes, "K0141403_1504", siblings.all, 
            chromlength)
ggsave(paste0("out/offsprings/", "K0141403_1504", "_with_offsprings.png"), 
       plot = gg, width = 8, height = 5)
gg <- familymapv2(chrombarcodes, "K0374502_1412", siblings.all, 
            chromlength)
ggsave(paste0("out/offsprings/", "K0374502_1412", "_with_offsprings.png"), 
       plot = gg, width = 8, height = 5)

##########################

#############Super family#############

superfamily.all <- c("K0374413_1611",
                     "K0060617_1512", 
                     "P0020206_1702",
                     "P0020206_1703",
                     "P0060623_1612",
                     "J0272804_1511",
                     "P0020215_1701")
for (superfam in superfamily.all){
  gg <- familymapv2(chrombarcodes, superfam, superfamily.all, 
                    chromlength)
  ggsave(paste0("out/siblings/distant/", superfam, "_with_distant-siblings.png"), 
         plot = gg, width = 8, height = 5)
}
##########################

#############IBD hotspots#############
conserved_regions <- read.csv("rawdata/Nwakanma-2013-ConservedRegions.csv")
colnames(conserved_regions) <- c("chr", "start", "end", "size", "SNPs", "Genes", "GeneNames")
conserved_regions$start <- conserved_regions$start*1000
conserved_regions$end <- conserved_regions$end*1000

edges.notsimilar <- edges[edges$fract_sites_IBD < 0.5,]
chrombarcodes.notsimilar <- chrombarcodes[paste0(chrombarcodes$ID1, chrombarcodes$ID2) %in% paste0(edges.notsimilar$ID1, edges.notsimilar$ID2),]
chrombarcodes.notsimilar.same <- chrombarcodes.notsimilar[chrombarcodes.notsimilar$different == 0,]
subj <-
  with(chrombarcodes.notsimilar.same, GenomicRanges::GRanges(chr, IRanges(start, end)))
chrombarcodes.notsimilar.same$coverage <- GenomicRanges::countOverlaps(subj, subj, type = "within")
chrombarcodes.notsimilar.same.unique <- chrombarcodes.notsimilar.same[!duplicated(paste(chrombarcodes.notsimilar.same$chr, chrombarcodes.notsimilar.same$start, 
                                                                                       chrombarcodes.notsimilar.same$end)),]
quantile(chrombarcodes.notsimilar.same$coverage, probs = c(0,0.05,0.1,0.25,0.5,0.75,0.9,0.95,1))

colorpalette <- brewer.pal(n = 11, name = "RdYlGn")
shareranges <- round(exp(seq(log(10),log(2190), length.out = 12))/10)*10
shareranges.df <- data.frame(shareranges[-length(shareranges)], shareranges[-1])
colnames(shareranges.df) <- c("from", "to")
shareranges.df$location <- seq(2,3,length.out = 11)
shareranges.df$color <- colorpalette

gg <- ggplot()+
  geom_rect(data=chromlength,
            aes(xmin = 0, xmax=length, 
                ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
            color=NA, fill="grey80")+
  geom_text(data=chromlength,
            aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                label=chr), 
            color="grey50", size=3)

for (i in c(1:11)){
  sharerange <- c(shareranges[i], shareranges[i+1])
  currcolor <- colorpalette[i]
  gg <- gg+
  geom_rect(data=chrombarcodes.notsimilar.same.unique[chrombarcodes.notsimilar.same.unique$coverage >= sharerange[1] & chrombarcodes.notsimilar.same.unique$coverage < sharerange[2],],
            aes(xmin=start, xmax=end, 
                ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
            fill = currcolor, size=0.1)
}

gg <- gg+
  geom_rect(data=conserved_regions,
             aes(xmin = start, xmax = end,
                 ymin=chr*2.75-0.6, ymax=chr*2.75+0.6),
            color = "blue", fill = NA, linetype = 5, size = 0.4)+
  geom_rect(data = data.frame(1),
            xmin = 1.35*10^6, xmax = 3.15*10^6,
            ymin = 3*2.75-1.5, ymax = 3*2.75+0.75, color = "black", 
            fill = "white", size = 0.5 )+
  geom_text(data = data.frame(1),
            x = 1.675*10^6, 
            y = 3*2.75, color = "black", 
            label = "Maximal pairwise identity occurence", size = 3.5)+
  geom_text(data = shareranges.df,
               aes(x = location*10^6,
                   y = 3*2.75-0.95, label = from), 
               color = "black", size = 3.5)+
  geom_text(data = shareranges.df,
               aes(x = location*10^6+(1*10^6/10),
                   y = 3*2.75-0.95, label = to), 
               color = "black", size = 3.5)+
  geom_segment(data = shareranges.df,
               aes(x = location*10^6+(1*10^6/10), xend = location*10^6+(1*10^6/10),
                   y = 3*2.75-0.75, yend = 3*2.75-0.5), 
               color = "black", size = 0.5)+
  geom_segment(data = shareranges.df,
               aes(x = location*10^6, xend = location*10^6,
                   y = 3*2.75-0.75, yend = 3*2.75-0.5), 
               color = "black", size = 0.5)+
  geom_rect(data=shareranges.df,
            aes(xmin=location*10^6, xmax=location*10^6+(1*10^6/10), 
                ymin=3*2.75-0.5, ymax=3*2.75+0.5, fill = color),
            size=0.1, color = NA)+
  scale_fill_manual(breaks = shareranges.df$color,
                    values = shareranges.df$color)+
  scale_x_continuous(breaks = seq(0,3.25*10^6,0.1*10^6), 
                     labels = seq(0,3.25*10^6,0.1*10^6)/10^6,
                     minor_breaks = seq(0,3.25*10^6,0.025*10^6),
                     expand = c(0.01,0.05))+
  scale_y_discrete(expand = c(0.01,0.01))+
  theme(plot.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey60", linetype = 3, size = 0.5),
        panel.grid.minor.x = element_line(color = "grey80", linetype = 3, size = 0.5),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "grey60", size = 0.5))
gg

ggsave(paste0("out/conserved-regions.png"), 
       plot = gg, width = 16, height = 10, dpi = 200)
##########################


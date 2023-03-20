#07/10/2021
#!/usr/bin/env Rscript

#This script will estimate the minimal allele frequency of the majority
#base called necessary to discard the minority base called as noise

#It will comparge WGs barcodes to genotyped barcodes and count the number of
#matches (mixed to mixed) and mismacthes (mixed to not mixed) resulting from
#the WGS barcodes generated using different values of cutoffs

###########################SCRIPT###########################
nnumbers <- data.frame(ID = NA, alignment = NA)
#All these cutoffs will be tested
minpropswgs <- seq(0.5,1,0.01)
for (minprop in minpropswgs){
  source("WGS-barcodes.R")
  source("merge-SNP-WGS.R")
  nnumbers <- merge(nnumbers, nnumber, all = T)
}


nnumbers <- nnumbers[!is.na(nnumbers$ID),]

###########################
###########################2. DETERMINE THE CUTOFF###########################

#We chose the cutoff maximizing N matches and minimizing N mismacthes
library(dplyr)
nnumbers.melted <- melt(nnumbers)
nnumbers.summary <- nnumbers.melted %>%
  group_by(alignment, variable) %>%
  summarise(Mean = mean(value))
#nn = N-N
#nl = N-ATGC
gg <- ggplot(data = nnumbers.summary[nnumbers.summary$alignment%in%c("nn", "nl"),], 
             aes(x = 1-as.numeric(sub("ratio", "", variable)), y = Mean, color = alignment))+
  geom_path(aes(group = alignment), size = 1)+
  geom_point(data = nnumbers.summary[nnumbers.summary$alignment%in%c("nn", "nl") & 
                                       nnumbers.summary$variable=="ratio0.8",],
             size = 3)+
  xlab("Threshold of within-sample MAF to call a genetic locus mixed")+
  ylab("Molecular vs genetic barcode mixed locus comparison")+
  scale_color_manual(name = "", labels = c("Mean number of mixed locus matching", "Mean number of mixed locus not matching"),
                     values = c("nn" = "cyan3", "nl" = "yellow2"))+
  scale_y_continuous(n.breaks = 15, expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), n.breaks = 15)+
  theme(axis.text.x = element_text(angle=90),
        legend.position = "top",
        text = element_text(size = 9))
gg
#This corresponds to 0.8 of minimal majority allele called
ggsave("out/wgss-barcode-cutoff.png", width = 8, height = 3)

rm(minprop)
######################################################
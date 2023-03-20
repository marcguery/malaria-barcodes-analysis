#####################READ DATA#################################
barcode.nodes <- read.table("../network/out/nodes.tsv", header = TRUE)
barcode.edges <- read.table("../network/out/edges.tsv", header = TRUE)

wgs.nodes <- read.table("../triades/out/nodes.tsv", header = TRUE)
wgs.edges <- read.table("../triades/out/edges.tsv", header = TRUE)
######################################################

#####################MERGE#################################
colnames(wgs.edges)[which(colnames(wgs.edges)=="fract_sites_IBD")] <- "fract_sites_IBD.wgs"
colnames(barcode.edges)[which(colnames(barcode.edges)=="fract_sites_IBD")] <- "fract_sites_IBD.barcode"

edges <- merge(barcode.edges, wgs.edges,
               suffixes = c(".barcode", ".wgs"))
nodes <- merge(barcode.nodes, wgs.nodes)
######################################################

#####################PLOT#################################
rsq <- function(x, y) summary(lm(y~x))$r.squared
lm05 <- lm(edges$fract_sites_IBD.barcode[edges$fract_sites_IBD.barcode > 0.5 & edges$fract_sites_IBD.wgs > 0.5] ~ 
     edges$fract_sites_IBD.wgs[edges$fract_sites_IBD.barcode > 0.5 & edges$fract_sites_IBD.wgs > 0.5])

rsq05 <- rsq(edges$fract_sites_IBD.wgs[edges$fract_sites_IBD.barcode > 0.5 & edges$fract_sites_IBD.wgs > 0.5],
    edges$fract_sites_IBD.barcode[edges$fract_sites_IBD.barcode > 0.5 & edges$fract_sites_IBD.wgs > 0.5])

accuracy <- data.frame(cutoff = c(),
                       sensitivity = c(), specificity = c(),
                       ppv = c(), npv = c(),
                       pos = c(), neg = c())
nrow(edges[edges$fract_sites_IBD.barcode>=0 & edges$fract_sites_IBD.wgs >= 0,])

ctoff <- 0.5
  tp <- length(which(edges$fract_sites_IBD.barcode >= ctoff & edges$fract_sites_IBD.wgs >= ctoff))
  fp <- length(which(edges$fract_sites_IBD.barcode >= ctoff & edges$fract_sites_IBD.wgs >= 0 & edges$fract_sites_IBD.wgs < ctoff))
  fn <- length(which(edges$fract_sites_IBD.barcode >= 0 & edges$fract_sites_IBD.barcode < ctoff & edges$fract_sites_IBD.wgs >= ctoff))
  tn <- length(which(edges$fract_sites_IBD.barcode>=0 & edges$fract_sites_IBD.barcode < ctoff & edges$fract_sites_IBD.wgs >= 0 & edges$fract_sites_IBD.wgs < ctoff))
  se <- tp/(tp+fn)
  sp <- tn/(fp+tn)
  ppv <- tp/(tp+fp)
  npv <- tn/(fn+tn)
  pos <- tp + fp
  neg <- tn + fn
  accuracy <- rbind(accuracy, data.frame(cutoff = ctoff,
                                         sensitivity = se,
                                         specificity = sp,
                                         ppv = ppv,
                                         npv = npv,
                                         pos = pos,
                                         neg = neg))

edges$accuracy <- NA
edges$accuracy[edges$fract_sites_IBD.barcode >= ctoff & edges$fract_sites_IBD.wgs >= ctoff] <- "tp"
edges$accuracy[edges$fract_sites_IBD.barcode >= ctoff & edges$fract_sites_IBD.wgs >= 0 & edges$fract_sites_IBD.wgs < ctoff] <- "fp"
edges$accuracy[edges$fract_sites_IBD.barcode >= 0 & edges$fract_sites_IBD.barcode < ctoff & edges$fract_sites_IBD.wgs >= ctoff] <- "fn"
edges$accuracy[edges$fract_sites_IBD.barcode>=0 & edges$fract_sites_IBD.barcode < ctoff & edges$fract_sites_IBD.wgs >= 0 & edges$fract_sites_IBD.wgs < ctoff] <- "tn"

gg1 <- ggplot(edges[edges$fract_sites_IBD.barcode>=0 & edges$fract_sites_IBD.wgs >= 0,])+
  geom_point(aes(x = fract_sites_IBD.wgs, y = fract_sites_IBD.barcode,
                  color = accuracy),
              show.legend = T)+
  scale_color_manual(name = "Classification",
                        values = c("tp" = "green4", "tn" = "red4",
                                  "fp" = "green2", "fn" = "red"),
                     labels = c("Truly related", "Truly unrelated",
                                "Falsely related", "Falsely unrelated"))+
  geom_segment(data=data.frame(1),
               x = 0, y = 0,
               xend = Inf, yend = Inf, linetype = 2, alpha = 0.8, size = 0.75)+
  xlab("Filtered IBD from WGS SNPs")+
  ylab("Filtered IBD from barcodes")+
  scale_x_continuous(limits = c(0,1), expand = c(0.01,0))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01,0))
  
gg1
ggsave("out/ibd-correlation.png",
       width = 8, height = 6)
######################################################

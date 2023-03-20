######################READING######################
rawdir <- "../read/rawdata/"
mappname <- "WGS-19262snps_mapping.csv"
vcfname <- "WGS-19262snps.vcf"
adname <- "WGS-19262snps_allele-depth.tsv"
samplestoremove <- read.table("../read/rawdata/samples-to-remove.txt", sep ="\t", header = TRUE)

mapping <- read.csv(paste0(rawdir, "../out/", mappname))

#gzip broken pipe expected because only the header is fetched
vcf_header <- system(str_glue('zcat "{rawdir}{vcfname}".zip | grep -m1 "#CHROM"'), ignore.stderr = T, intern = T)
vcf_header <- strsplit(vcf_header, split = "\t")[[1]]
snps <- read.delim(pipe(paste("zcat", paste0(rawdir, vcfname, ".zip"))), 
                           comment.char = "#", header = FALSE)
colnames(snps) <- make.names(vcf_header)
snps <- snps[,c(1:9,9+which(!colnames(snps[,10:ncol(snps)])%in%sub(".WGS", "", samplestoremove$ID)))]
snps.ad <- read.delim(pipe(paste("zcat", paste0(rawdir, adname, ".zip"))), 
                      header = TRUE)
snps.ad <- snps.ad[,c(1:2,2+which(!colnames(snps.ad[,3:ncol(snps.ad)])%in%sub(".WGS", "", samplestoremove$ID)))]
############################################

######################FORMATTING######################
snps <- snps[,c(1:5)]
colnames(snps)[1] <- "CHROM"

snps.merged <- merge(snps, snps.ad, by = c("CHROM", "POS"))
snps.merged <- snps.merged[snps.merged$CHROM!="Pf3D7_API_v3",]

snps.ref <- as.data.frame(snps.merged[,-c(1:5)])
snps.ref <- as.data.frame(apply(snps.ref, 2, function(x){as.integer(sub(",\\S+", "", x))}),
                          row.names=paste(snps.merged[,1], snps.merged[,2], sep=":"))

snps.alt <- as.data.frame(snps.merged[,-c(1:5)], 
                          row.names = paste(snps.merged[,1], snps.merged[,2], sep=":"))
snps.alt <- as.data.frame(apply(snps.alt, 2, function(x){as.integer(sub("\\S+,", "", x))}),
                          row.names=paste(snps.merged[,1], snps.merged[,2], sep=":"))

############################################
#################CALCULATE HETERO#################

lowdepth <- 5
snpsnum <- nrow(snps.merged)
minsnpnum <- 4000

snps.sum <- snps.ref+snps.alt
snps.ref <- snps.ref[,which(apply(snps.sum, 2, function(x){
  length(which(x>lowdepth))})>minsnpnum)]
snps.alt <- snps.alt[,which(apply(snps.sum, 2, function(x){
  length(which(x>lowdepth))})>minsnpnum)]

snps.sum <- snps.ref+snps.alt
snps.ref[snps.sum<=lowdepth] <- 0
snps.alt[snps.sum<=lowdepth] <- 0
snps.ref[snps.sum==0] <- NA
snps.alt[snps.sum==0] <- NA

p.sample<-snps.ref/(snps.ref + snps.alt)
q.sample<-snps.alt/(snps.ref + snps.alt)
h.sample<-1-(p.sample^2 + q.sample^2)  # h = 1- (p^2 / q^2)      " h " represents the heterozygosity for a particular SNP in a particular sample

# Calculate heterozygosity for each SNP across all samples as 1-(F1^2 + F2^2)
snps.ref.sums<-apply(snps.ref,1,sum,na.rm=T)
snps.alt.sums<-apply(snps.alt,1,sum,na.rm=T)
P.pop<-snps.ref.sums/(snps.ref.sums + snps.alt.sums)    #
Q.pop<-snps.alt.sums/(snps.ref.sums + snps.alt.sums)   # a.k.a, 1-F1   #  F1 <- 1-F2
H.pop<-1 - (P.pop^2 + Q.pop^2)

Fboth<-data.frame(P.pop, Q.pop)
MAF<-apply(Fboth,1,min)
write.table(MAF, "out/maf.tsv", row.names = T, col.names = F, quote = F, sep ="\t")
# We want to bin the SNPs into 10 bins in order of MAF
bins=list()
for (num in seq(0,0.45, 0.05)){
  bins[[paste("l", as.character(num), "h", as.character(num+0.05), sep="")]]<-which(MAF>num & MAF<=(num+0.05))
}

h.10 <- bins
h.10 <- lapply(h.10, function(x){apply(h.sample[x,,drop=F], 2, mean, na.rm=T)})

H.10 <- bins
H.10 <- lapply(H.10, function(x){mean(H.pop[x])})

h.10 <- t(as.data.frame(h.10))
H.10 <- unlist(H.10)

h.10[is.nan(h.10)] <- 0
H.10[is.nan(H.10)] <- 0

b.4<-rep(0,ncol(h.10)) # This will hold our beta's

for(i in 1:ncol(h.10)){
  b.4[i]<-lm(h.10[,i] ~ H.10 -1)$coefficients[1]
}

# Calculate inbreeding value for each sample as in the methods
Fws<- 1-b.4
names(Fws) <- colnames(h.sample)
Fws <- data.frame(Fws)
snps.sum <- snps.ref+snps.alt
snpsUsed <- data.frame(apply(snps.sum, 2, function(x){
  length(which(x>lowdepth))}))
colnames(snpsUsed) <- "SNPs"
Fws <- merge(Fws, snpsUsed, by="row.names")

mapping$sample <- make.names(mapping$sample)
mapping$sample <- sub(".WGS$", "", mapping$sample)
mapping <- mapping[!is.na(mapping$SampleID_corrected),]
#Removal of one bad format sample (PA0761.CW) and one 3D7 related sample (SPT21929)
Fws <- merge(Fws, mapping, by.x = "Row.names", by.y = "sample")
Fws <- Fws[,c("SampleID_corrected", "Fws", "SNPs")]
colnames(Fws) <- c("Sample", "Fws", "SNPs")
write.csv(Fws, paste("out/fws", snpsnum, lowdepth, minsnpnum, "persample.csv", sep="-"), quote = F, row.names = F)

##################################

#################PLOT#################
### To make plot as in the Manske2012 Nature paper

heteromat <- as.data.frame(cbind(H.10, h.10))
heteromat <- heteromat[,c(1,sample(2:ncol(heteromat),20))]
heteromat <- melt(heteromat, id.vars = c(1))
heteromat <- merge(heteromat, mapping, by.x = "variable", by.y = "sample")
heteromat <-  heteromat[,c("SampleID_corrected", "H.10", "value")]
colnames(heteromat) <- c("Sample", "PopH", "SampH")

gg <- ggplot(heteromat)+
  geom_smooth(aes(x=PopH, y=SampH, col=factor(Sample)),
              method='lm', formula= y~x, se = F)+
  geom_point(aes(x=PopH, y=SampH, col=factor(Sample)))+
  xlab("Population heterozygosity")+
  ylab("Sample Heterozygosity")
gg

ggsave("out/Fws-popH.png", width = 8, units = "in")

gg <- ggplot(Fws, aes(x = factor(-0.05+floor((Fws+0.05)*20)/20), 
                      y = Fws,
                      color = factor(-0.05+floor((Fws+0.05)*20)/20),
                      fill = factor(-0.05+floor((Fws+0.05)*20)/20),
                      group = factor(-0.05+floor((Fws+0.05)*20)/20)))+
  geom_col(show.legend = F)+
  scale_y_continuous(breaks = seq(0,nrow(Fws),0.6*nrow(Fws)/12), 
                     limits = c(0,0.605*nrow(Fws)),
                     labels = function(x){100*round(x/nrow(Fws),2)})+
  scale_x_discrete(labels = function(x){paste0(as.numeric(x),"-",as.numeric(x)+0.05)})+
  ylab("Samples (%)")+
  xlab("Fws")+
  theme(panel.grid.major.x = element_blank(),
        text = element_text(size = 18))
gg

ggsave("out/Fws-distribution.png", width = 15, height = 6, units = "in")
##################################

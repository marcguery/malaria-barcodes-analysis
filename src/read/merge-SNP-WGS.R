###############SNP and WGS BARCODES###############
WGSbarcodes <- read.csv("out/WGS-barcodes.csv")
SNPbarcodes <- read.csv("out/SNP-barcodes.csv")
##############################
##################CLEANED DATA##################
barcodes.used.snp <- read.csv("../../101-genotyped-SNPS.csv", header = TRUE)
barcodes.used.snp$Study <- "SNP"
barcodes.used.wgs <- read.csv("../../101-WGS-SNPS.csv", header = TRUE)
barcodes.used.wgs$Study <- "WGS"
barcodes.used <- rbind(barcodes.used.snp, barcodes.used.wgs)
##############################

###############MERGING BARCODES###############
WGSavail <- WGSbarcodes[,colnames(WGSbarcodes)%in%colnames(SNPbarcodes)]
SNPavail <- SNPbarcodes[,colnames(SNPbarcodes)%in%colnames(WGSbarcodes)]

barcodes <- rbind(WGSavail, SNPavail)

#startPosindex and endPosIndex delimit the locations of the columns storing
#the bases called by SNP genotyping or WGS
startPosindex <- 2
endPosIndex <- ncol(barcodes)-2
barcodes$Barcode <- apply(barcodes[,c(startPosindex:endPosIndex)], MARGIN = 1,
                          FUN = function(x) {paste(x, sep="", collapse = "")})
##############################

################SNP and WGS MAPPING#################
WGSmapping <- read.csv("out/WGS-barcodes_mapping.csv")
colnames(WGSmapping) <- paste(colnames(WGSmapping), "WGS", sep=".")
SNPmapping <- read.csv("out/SNP-barcodes_mapping.csv")
colnames(SNPmapping) <- paste(colnames(SNPmapping), "SNP", sep=".")
################################

################MERGING MAPPING################
barcodes <- merge(barcodes, WGSmapping[,c(1,3)], 
            by.x = "Internal.Sample.ID", by.y = "sample.WGS", all.x = T)
barcodes$SampleID_corrected.WGS[barcodes$Study=="SNP"] <- NA
barcodes <- merge(barcodes, SNPmapping[,c(1,3)],
                  by.x = "Internal.Sample.ID", by.y = "Internal.Sample.ID.SNP", all.x = T)
barcodes$Full_ID.SNP[barcodes$Study=="WGS"] <- NA

barcodes$ID <- as.vector(barcodes$Full_ID.SNP)
barcodes$ID[is.na(barcodes$ID)] <- as.vector(barcodes$SampleID_corrected.WGS[is.na(barcodes$ID)])

which(is.na(barcodes$ID)) #REF and ALT should not have an ID, 
                          #as well as bad formatted samples from SNP dataset
barcodes.refalt <- barcodes[which(barcodes$Internal.Sample.ID%in%c("REF.WGS", "ALT.WGS")),]
barcodes.refalt$ID[which(barcodes.refalt$Internal.Sample.ID=="REF.WGS")] <- "REF.WGS"
barcodes.refalt$ID[which(barcodes.refalt$Internal.Sample.ID=="ALT.WGS")] <- "ALT.WGS"
barcodes.refalt$Note <- "Artificial"
#We remove these samples
barcodes.na <- data.frame(barcodes[which(is.na(barcodes$ID)),"Internal.Sample.ID"])
barcodes.na$why <- "Bad format"
colnames(barcodes.na)[1] <- "ID"
#write.table(barcodes.na, "rawdata/samples-to-remove.txt", append = TRUE, quote = F, row.names = F, col.names = T)
barcodes <- barcodes[-which(is.na(barcodes$ID)),]
################################

################CLEANING################
#Removing worse of duplicates
remove_duplicate <- function(barcodeID, study){
  dupbarcode <- barcodes[barcodes$ID==barcodeID & !is.na(barcodes$ID) & barcodes$Study==study,]
  if (nrow(dupbarcode) == 1) {
    return()
  }
  scores <- sapply(dupbarcode[,"Barcode"],
        FUN = function(x) {nchar(gsub(pattern = "[X|-]", replacement = "", x = x))})
  return(dupbarcode[-which.max(scores),"Internal.Sample.ID"])
}

studies <- c("SNP", "WGS")
barcodes.dup <- barcodes[1,]
barcodes.dup <- apply(barcodes.dup, c(1,2), function(x){x <- NA})
for (stud in studies){
  #From one of the raw data file
  idToRm <- unlist(sapply(names(which(table(barcodes$ID)>1)), remove_duplicate, study=stud))
  #Number of duplicated samples
  length(idToRm)
  #Saving duplicates
  barcodes.dup <- rbind(barcodes.dup, barcodes[barcodes$Internal.Sample.ID%in%idToRm,])
  #Removing duplicates
  barcodes <- barcodes[!barcodes$Internal.Sample.ID%in%idToRm,]
}
barcodes.dup <- data.frame(barcodes.dup[!is.na(barcodes.dup$Internal.Sample.ID),"Internal.Sample.ID"])
barcodes.dup$why <- "Duplicated"
#write.table(barcodes.dup, "rawdata/samples-to-remove.txt", append = TRUE, col.names = F, row.names = F, quote = F)
barcodes.dedupl <- barcodes
colorder.dedupl <- colnames(barcodes.dedupl)
barcodes.dedupl <- merge(barcodes.dedupl, barcodes.used[,c("ID", "Study")],
                                all.y = TRUE)
barcodes.dedupl <- barcodes.dedupl[,colorder.dedupl]
################################

##########CORRELATE WGS and SNP barcodes##########

#Attribute a score of alignment between SNP and WGS barcodes from the same blood sample
positionqual <- function(pairIBD, bcodes){
  snpbcode <- bcodes[bcodes$Study != "WGS" & !is.na(bcodes$Full_ID.SNP) & bcodes$Full_ID.SNP==pairIBD,]
  wgsbcode <- bcodes[bcodes$Study == "WGS" & !is.na(bcodes$SampleID_corrected.WGS) & bcodes$SampleID_corrected.WGS==pairIBD,]
  bcode <- rbind(snpbcode, wgsbcode)
  bcode <- t(bcode)
  colnames(bcode) <- bcode[nrow(bcode),]
  bcode <- bcode[c(startPosindex:endPosIndex),]
  bcode[bcode=="A"] <- 2
  bcode[bcode=="T"] <- 4
  bcode[bcode=="C"] <- 6
  bcode[bcode=="G"] <- 8
  bcode[bcode=="N"] <- 0
  bcode[bcode=="X"] <- "NA"
  bcode <- cbind(bcode, sub("['.']\\S+", "", row.names(bcode)))
  bcode <- cbind(bcode, sub("\\S+['.']", "", row.names(bcode)))
  colnames(bcode)[(ncol(bcode)-1):ncol(bcode)] <- c("chr", "pos")
  bcode <- as.data.frame(bcode)
  bcode$chr <- sub(pattern = "^\\S{5}_(0)?", replacement="", bcode$chr)
  bcode$chr <- sub(pattern = "_\\S+", replacement="", bcode$chr)
  scores <- as.integer(bcode[,1])/as.integer(bcode[,2])
  scores[scores==Inf|scores==0] <- -2 #N-ATGC : e. g. 2/0 or 0/2 
  scores[scores>0 & scores!=1 & scores<Inf] <- -1 #Mismatch : e. g. 6/2, 2/6
  scores[is.nan(scores)] <- 2 #N-N : 0/0=nan
  scores[is.na(scores)] <- 0 #X-ATGCNX : e.g. 2/NA, NA/0 = NA
  names(scores) <- row.names(bcode)
  return(scores)
}

#Samples having a SNPs and WGS barcode
dups <- names(which(table(barcodes.dedupl$ID)==2))

##########COMPARING UNFILTERED BARCODES##########
#Code for each site comparison
#This will trigger a warning as there are expected NAs in the scores
bcodescores.dedupl <- sapply(dups,
                      FUN = positionqual, barcodes.dedupl)

bcodescores.dedupl.melted <- melt(bcodescores.dedupl)

#Summary of the comparison of each SNP
bcodescores.dedupl.summary <- as.data.frame(
  matrix(c(apply(bcodescores.dedupl, 1,
                 FUN = function(x) {length(which(x==1))/length(which(x!=0))}),
           apply(bcodescores.dedupl, 1,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.dedupl, 1,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.dedupl, 1,
                 FUN = function(x) {length(which(x==-1))}),
           apply(bcodescores.dedupl, 1,
                 FUN = function(x) {length(which(x==0))}),
           apply(bcodescores.dedupl, 1,
                 FUN = function(x) {length(which(x==2))}),
           apply(bcodescores.dedupl, 1,
                 FUN = function(x) {length(which(x==-2))})
  ),
  ncol = 7
  ),
  row.names = row.names(bcodescores.dedupl)
)
colnames(bcodescores.dedupl.summary) <- c("meanComp", "sumSame", "same", "different", "unknown", "nn", "nl")

#Summary of the comparison of each SNP in each barcode
bcodescores.dedupl.summaryperbarcode <- as.data.frame(
  matrix(c(apply(bcodescores.dedupl, 2,
                 FUN = function(x) {length(which(x==1))/length(which(x!=0))}),
           apply(bcodescores.dedupl, 2,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.dedupl, 2,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.dedupl, 2,
                 FUN = function(x) {length(which(x==-1))}),
           apply(bcodescores.dedupl, 2,
                 FUN = function(x) {length(which(x==0))}),
           apply(bcodescores.dedupl, 2,
                 FUN = function(x) {length(which(x==2))}),
           apply(bcodescores.dedupl, 2,
                 FUN = function(x) {length(which(x==-2))})
  ),
  ncol = 7
  ),
  row.names = row.names(bcodescores.dedupl)
)
colnames(bcodescores.dedupl.summaryperbarcode) <- c("meanComp", "sumSame", "same", "different", "unknown", "nn", "nl")

#############RESETTING BAD SNP OF SNP > 1605################

#PLOT ORDER based on SNP quality
#Order of SNPs for plotting: 
##least number of individuals with different SNPs --> most
##least number of individuals with unkowns SNPs --> most
orderedsnps <- rownames(bcodescores.dedupl.summary[with(bcodescores.dedupl.summary, order(-different, -unknown)),])
#Order of barcode comparison for plotting
##past --> present
##least number of different SNPs --> most
##least number of unkowns SNPs --> most
dateperbarcode <- sub("\\S+_", "", colnames(bcodescores.dedupl))
datetile <- sub("\\S+_", "", bcodescores.dedupl.melted$Var2)
bcodescores.dedupl.summaryperbarcode$date <- as.numeric(dateperbarcode)
bcodescores.dedupl.melted$date <- as.numeric(datetile)
orderedbarcodes <- colnames(bcodescores.dedupl[,with(bcodescores.dedupl.summaryperbarcode, order(-date, -different, -unknown))])

#SNP LOCATION FOR PLOTTING
#Length of chromosomes
chlen <- read.table("rawdata/chrom-sizes.tsv", header = F)
colnames(chlen) <- c("ID", "length")
chlen$chr <- sub("Pf3D7_(0|)", "", chlen$ID)
chlen$chr <- sub("_v3", "", chlen$chr)

#All SNPs and all available SNPs (after merging of WGS and genotyped SNPs)
allSNPs <- read.table("rawdata/SNP-positions.tsv", header = F)
if (file.exists("rawdata/SNP-available.tsv")){
  availSNPs <- read.table("rawdata/SNP-available.tsv", header = F)
}else{
  availSNPs <- read.table("rawdata/SNP-positions.tsv", header = F)
}
colnames(allSNPs) <- c("chr", "pos")
colnames(availSNPs) <- c("chr", "pos")
availSNPs$avail <- "yes"
snps <- merge(availSNPs, allSNPs, by=c("chr", "pos"), all=T)
snps$avail[is.na(snps$avail)] <- "no"
snps$chr <- sub(pattern = "^\\S{5}_(0)?", replacement="", snps$chr)
snps$chr <- sub(pattern = "_\\S+", replacement="", snps$chr)

bcodescores.dedupl.summary$chr <- sub(pattern = "^\\S{5}_(0)?", replacement="", row.names(bcodescores.dedupl.summary))
bcodescores.dedupl.summary$pos <- sub(pattern = "^\\S+['.']", replacement="", bcodescores.dedupl.summary$chr)
bcodescores.dedupl.summary$chr <- sub(pattern = "_\\S+", replacement="", bcodescores.dedupl.summary$chr)
snps.dedupl <- merge(snps, bcodescores.dedupl.summary, all.x=T)

#Statistics of SNPs on each chromosome
snpusedperchs <- table(availSNPs$chr)
snpusedperchs <- as.data.frame(snpusedperchs)
colnames(snpusedperchs) <- c("chr", "Usedsnp")
snpperchs <- table(allSNPs$chr)
snpperchs <- as.data.frame(snpperchs)
colnames(snpperchs) <- c("chr", "Tsnp")
snpsnumber <- merge(snpperchs, snpusedperchs)
snpsnumber$chr <- sub(pattern = "^\\S{5}_(0)?", replacement="", snpsnumber$chr)
snpsnumber$chr <- sub(pattern = "_\\S+", replacement="", snpsnumber$chr)

chs <- merge(snpsnumber, chlen)

#Removal of 21 worse SNPs after May 2016
worsesnps.heatmap.merged <- orderedsnps[1:21]
worsesnps.heatmap <- data.frame(matrix(nrow=21))
worsesnps.heatmap[,1] <- sub("\\.\\S+", "", worsesnps.heatmap.merged)
worsesnps.heatmap$pos <- sub("\\S+\\.", "", worsesnps.heatmap.merged)
colnames(worsesnps.heatmap) <- c("chr", "pos")
worsesnps.heatmap$chr <- sub("_v3", "", worsesnps.heatmap$chr)
worsesnps.heatmap$chr <- sub("^\\S+_", "", worsesnps.heatmap$chr)
worsesnps.heatmap$chr <- as.integer(worsesnps.heatmap$chr)
worsesnps.comp <- snps.dedupl[order(snps.dedupl$meanComp)[1:21], c(1,2)]
worsesnps <- merge(worsesnps.comp, worsesnps.heatmap, all = T)
nrow(worsesnps)==length(worsesnps.heatmap.merged)
#              ^ Should be True
#SNPs that are not showing the same information between WGS and genotyped barcodes
worsesnps.heatmap.merged

#Removal of the discordant SNPs
barcodesdate <- sub("\\S+_", "", barcodes$ID)

barcodes[as.integer(barcodesdate) > 1605 & barcodes$Study!="WGS", which(colnames(barcodes)%in%worsesnps.heatmap.merged)] <- "X"
barcodes.goodsnps <- barcodes

colorder.goodsnps <- colnames(barcodes.goodsnps)
barcodes.goodsnps <- merge(barcodes.goodsnps, barcodes.used[,c("ID", "Study")],
                       all.y = TRUE)
barcodes.goodsnps <- barcodes.goodsnps[,colorder.goodsnps]

#############COMPARING BARCODES WITHOUT BAD SNPS#############
#ESTIMATE THE CUTOFF FOR N in WGS barcodes
#Barcodes without bad SNPs
#This will trigger a warning as there are expected NAs in the scores
bcodescores.goodsnps <- sapply(dups,
                             FUN = positionqual, barcodes.goodsnps)

bcodescores.goodsnps.melted <- melt(bcodescores.goodsnps)

#Summary of the comparison of each SNP
bcodescores.goodsnps.summary <- as.data.frame(
  matrix(c(apply(bcodescores.goodsnps, 1,
                 FUN = function(x) {length(which(x==1))/length(which(x!=0))}),
           apply(bcodescores.goodsnps, 1,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.goodsnps, 1,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.goodsnps, 1,
                 FUN = function(x) {length(which(x==-1))}),
           apply(bcodescores.goodsnps, 1,
                 FUN = function(x) {length(which(x==0))}),
           apply(bcodescores.goodsnps, 1,
                 FUN = function(x) {length(which(x==2))}),
           apply(bcodescores.goodsnps, 1,
                 FUN = function(x) {length(which(x==-2))})
  ),
  ncol = 7
  ),
  row.names = row.names(bcodescores.goodsnps)
)
colnames(bcodescores.goodsnps.summary) <- c("meanComp", "sumSame", "same", "different", "unknown", "nn", "nl")

#Summary of the comparison of each SNP in each barcode
bcodescores.goodsnps.summaryperbarcode <- as.data.frame(
  matrix(c(apply(bcodescores.goodsnps, 2,
                 FUN = function(x) {length(which(x==1))/length(which(x!=0))}),
           apply(bcodescores.goodsnps, 2,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.goodsnps, 2,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.goodsnps, 2,
                 FUN = function(x) {length(which(x==-1))}),
           apply(bcodescores.goodsnps, 2,
                 FUN = function(x) {length(which(x==0))}),
           apply(bcodescores.goodsnps, 2,
                 FUN = function(x) {length(which(x==2))}),
           apply(bcodescores.goodsnps, 2,
                 FUN = function(x) {length(which(x==-2))})
  ),
  ncol = 7
  ),
  row.names = colnames(bcodescores.goodsnps)
)
colnames(bcodescores.goodsnps.summaryperbarcode) <- c("meanComp", "sumSame", "same", "different", "unknown", "nn", "nl")

#Adding dates
dateperbarcode <- sub("\\S+_", "", colnames(bcodescores.goodsnps))
datetile <- sub("\\S+_", "", bcodescores.goodsnps.melted$Var2)
bcodescores.goodsnps.summaryperbarcode$date <- as.numeric(dateperbarcode)
bcodescores.goodsnps.melted$date <- as.numeric(datetile)
#Snp discordance between genotyped SNP and WGS SNP
bcodescores.goodsnps.summary$chr <- sub(pattern = "^\\S{5}_(0)?", replacement="", row.names(bcodescores.goodsnps.summary))
bcodescores.goodsnps.summary$pos <- sub(pattern = "^\\S+['.']", replacement="", bcodescores.goodsnps.summary$chr)
bcodescores.goodsnps.summary$chr <- sub(pattern = "_\\S+", replacement="", bcodescores.goodsnps.summary$chr)
snps.goodsnps <- merge(snps, bcodescores.goodsnps.summary, all.x=T)

####################################

#############DIFFERENCES IN NUMBER OF MIXED CALLS BETWEEN WGS AND BARCODES##############
#This dataframe can be used to determine the appropriate cutoff
#of the minimal allele frequency of the majority base called
#to discard the minority base called
#(i. e., to add mixed calls inWGS barcodes)
bcodescores.goodsnps.summaryperbarcode$ID <- row.names(bcodescores.goodsnps.summaryperbarcode)
nnumber <- melt(bcodescores.goodsnps.summaryperbarcode[,c("ID", "same", "different", "unknown", "nn", "nl")])
if (!exists("minprop")){
  minprop <- ""
}
colnames(nnumber) <- c("ID", "alignment", paste0("ratio", minprop))
############################

##############MAKING CONSENSUS BARCODES###############
#Consensus barcodes are the same as original barcodes with:
#Barcodes from WGS replace all unknown calls of genotyped barcodes
#Sites with a discrepancy between SNP and WGS calls are set to unkown call
make_consensus <- function(barcodeID){
  ngbarcode <- barcodes[barcodes$ID==barcodeID & !is.na(barcodes$ID),]
  nbarcode <- ngbarcode[ngbarcode$Study!="WGS",c(startPosindex:endPosIndex)]
  gbarcode <- ngbarcode[ngbarcode$Study=="WGS",c(startPosindex:endPosIndex)]
  #Replace SNP Xs by WGS call
  cbarcode <- nbarcode
  cbarcode[which(cbarcode=="X")] <- gbarcode[which(cbarcode=="X")]
  cbarcode[which(cbarcode=="N")] <- gbarcode[which(cbarcode=="N")]
  cbarcode[which(cbarcode!="X" & cbarcode!="N" & gbarcode!="X" & gbarcode!="N" & cbarcode != gbarcode)] <- "X"
  res <- as.vector(t(cbarcode))
  names(res) <- colnames(cbarcode)
  return(res)
}

#Updating barcodes with consensus barcodes, removing duplicated samples from
#SNP genotyping and WGS
barcodes.dupl <- as.data.frame(t(sapply(names(which(table(barcodes$ID)==2)), make_consensus)))
barcodes[match(row.names(barcodes.dupl),barcodes$Full_ID.SNP),c(startPosindex:endPosIndex)] <- barcodes.dupl
barcodes$Note <- "Raw"
barcodes[match(row.names(barcodes.dupl),barcodes$Full_ID.SNP),"Note"] <- "Consensus"
barcodes$Barcode <- apply(barcodes[,c(startPosindex:endPosIndex)], MARGIN = 1,
                          FUN = function(x) {paste(x, sep="", collapse = "")})
barcodes.cons <- barcodes
colorder.cons <- colnames(barcodes.cons)
barcodes.cons <- merge(barcodes.cons, barcodes.used[,c("ID", "Study")],
                         all.y = TRUE)
barcodes.cons <- barcodes.cons[,colorder.cons]

#This will trigger a warning as there are expected NAs in the scores
bcodescores.cons <- sapply(dups,
                               FUN = positionqual, barcodes.cons)

bcodescores.cons.melted <- melt(bcodescores.cons)

#Summary of the comparison of each SNP
bcodescores.cons.summary <- as.data.frame(
  matrix(c(apply(bcodescores.cons, 1,
                 FUN = function(x) {length(which(x==1))/length(which(x!=0))}),
           apply(bcodescores.cons, 1,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.cons, 1,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.cons, 1,
                 FUN = function(x) {length(which(x==-1))}),
           apply(bcodescores.cons, 1,
                 FUN = function(x) {length(which(x==0))}),
           apply(bcodescores.cons, 1,
                 FUN = function(x) {length(which(x==2))}),
           apply(bcodescores.cons, 1,
                 FUN = function(x) {length(which(x==-2))})
  ),
  ncol = 7
  ),
  row.names = row.names(bcodescores.cons)
)
colnames(bcodescores.cons.summary) <- c("meanComp", "sumSame", "same", "different", "unknown", "nn", "nl")

#Summary of the comparison of each SNP in each barcode
bcodescores.cons.summaryperbarcode <- as.data.frame(
  matrix(c(apply(bcodescores.cons, 2,
                 FUN = function(x) {length(which(x==1))/length(which(x!=0))}),
           apply(bcodescores.cons, 2,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.cons, 2,
                 FUN = function(x) {length(which(x==1))}),
           apply(bcodescores.cons, 2,
                 FUN = function(x) {length(which(x==-1))}),
           apply(bcodescores.cons, 2,
                 FUN = function(x) {length(which(x==0))}),
           apply(bcodescores.cons, 2,
                 FUN = function(x) {length(which(x==2))}),
           apply(bcodescores.cons, 2,
                 FUN = function(x) {length(which(x==-2))})
  ),
  ncol = 7
  ),
  row.names = colnames(bcodescores.cons)
)
colnames(bcodescores.cons.summaryperbarcode) <- c("meanComp", "sumSame", "same", "different", "unknown", "nn", "nl")

#Adding dates
dateperbarcode <- sub("\\S+_", "", colnames(bcodescores.cons))
datetile <- sub("\\S+_", "", bcodescores.cons.melted$Var2)
bcodescores.cons.summaryperbarcode$date <- as.numeric(dateperbarcode)
bcodescores.cons.melted$date <- as.numeric(datetile)
#Snp discordance between genotyped SNP and WGS SNP
bcodescores.cons.summary$chr <- sub(pattern = "^\\S{5}_(0)?", replacement="", row.names(bcodescores.cons.summary))
bcodescores.cons.summary$pos <- sub(pattern = "^\\S+['.']", replacement="", bcodescores.cons.summary$chr)
bcodescores.cons.summary$chr <- sub(pattern = "_\\S+", replacement="", bcodescores.cons.summary$chr)
snps.cons <- merge(snps, bcodescores.cons.summary, all.x=T)


#Removing WGS barcodes having a new consensus
barcodes <- barcodes[-which(barcodes$Study=="WGS" & barcodes$ID%in%row.names(barcodes.dupl)),]

##############FILTERING LOW QUALITY BARCODES##############
#Removing X only barcodes
goodsnps <- sapply(X = barcodes$Barcode,
                   FUN = function(x) {nchar(gsub(pattern = "[X|-]", replacement = "", x = x))})
#Distribution of number of good quality SNPs in barcodes
quantile(goodsnps)
#Percentage of barcodes having no good SNPs at all
ecdf(goodsnps)(0)*100
#number of barcodes having no good SNPs at all
length(which(goodsnps==0))

#Removing samples having no good snps
barcodes <- barcodes[goodsnps > 0,]

barcodes.original <- barcodes

length(unique(barcodes$ID))
#This should be nrow
goodsnps <- sapply(X = barcodes$Barcode,
                   FUN = function(x) {nchar(gsub(pattern = "[X|-]", replacement = "", x = x))})
min(goodsnps) #This should be greater than 0
#We need to get rid of barcodes having too much unknown SNPs
ctoffsnps <- 21 #Minimal number of SNPs for a barcode to be kept
#Distribution of number of good quality SNPs in barcodes
quantile(goodsnps)
#Percentage of barcodes having no good SNPs at all
ecdf(goodsnps)(ctoffsnps-1)*100
#number of barcodes having no good SNPs at all
length(which(goodsnps<ctoffsnps))

barcodes <- barcodes[goodsnps >= ctoffsnps,]
barcodes.filtered <- barcodes
barcodes.refalt.filtered <- rbind(barcodes.filtered, barcodes.refalt)
############################

##############STATISTICS##############
#DEDUPL
meta.dedupl <- barcodes.dedupl[,c("Study", "ID")]
meta.dedupl$Study <- sub("\\d", "", meta.dedupl$Study)

meta.dedupl$Date <- sub("\\S+_", "", meta.dedupl$ID)
meta.dedupl$Date <- as.Date(paste0(meta.dedupl$Date, "01"), "%y%m%d")
meta.dedupl$Village <- toupper(substr(meta.dedupl$ID, 1, 1))
meta.dedupl <- meta.dedupl[meta.dedupl$Village%in%c("J", "K", "P", "N"),]

#GOODSNPS
meta.goodsnps <- barcodes.goodsnps[,c("Study", "ID")]
meta.goodsnps$Study <- sub("\\d", "", meta.goodsnps$Study)

meta.goodsnps$Date <- sub("\\S+_", "", meta.goodsnps$ID)
meta.goodsnps$Date <- as.Date(paste0(meta.goodsnps$Date, "01"), "%y%m%d")
meta.goodsnps$Village <- toupper(substr(meta.goodsnps$ID, 1, 1))
meta.goodsnps <- meta.goodsnps[meta.goodsnps$Village%in%c("J", "K", "P", "N"),]

#CONS
meta.cons <- barcodes.cons[,c("Study", "ID")]
meta.cons$Study <- sub("\\d", "", meta.cons$Study)

meta.cons$Date <- sub("\\S+_", "", meta.cons$ID)
meta.cons$Date <- as.Date(paste0(meta.cons$Date, "01"), "%y%m%d")
meta.cons$Village <- toupper(substr(meta.cons$ID, 1, 1))
meta.cons <- meta.cons[meta.cons$Village%in%c("J", "K", "P", "N"),]

#FILTERED
meta.filtered <- barcodes.filtered[,c("Study", "ID")]
meta.filtered$Study <- sub("\\d", "", meta.filtered$Study)

meta.filtered$Date <- sub("\\S+_", "", meta.filtered$ID)
meta.filtered$Date <- as.Date(paste0(meta.filtered$Date, "01"), "%y%m%d")
meta.filtered$Village <- toupper(substr(meta.filtered$ID, 1, 1))
meta.filtered <- meta.filtered[meta.filtered$Village%in%c("J", "K", "P", "N"),]

#SEASONS

#dry (january to july) wet (beginning august to end december)
wetseason <- c(8,12)
wetseason <- as.character(wetseason)

wetseason[nchar(wetseason)==1] <- paste0("0", wetseason[nchar(wetseason)==1])

years <- sort(as.numeric(unique(format(meta.dedupl$Date,"%Y"))))
wetseason <- as.Date(c(paste0(years[1], "-01-01"),
                       unlist(lapply(years, function(x){
                         paste(x, wetseason, "15", sep="-")
                       })),
                       paste0(years[length(years)], "-12-31")), "%Y-%m-%d")

seasons <- data.frame(wetseason[-length(wetseason)], wetseason[-1])
seasons$type <- rep(c("dry", "wet"), length.out=nrow(seasons))
colnames(seasons) <- c("date1", "date2", "type")
mindate <- min(meta.dedupl$Date)
#maxdate <- max(plotdataset$Date)
maxdate <- as.Date("2017-01-01")
#Time zone warning but it is ok
seasons <- seasons[seasons$date2>mindate & seasons$date1 < maxdate,]

yearspan <- data.frame(as.Date("2014-01-01")+months(c(0,12,24,36)))
colnames(yearspan) <- "year"
yearspan$manualposition <- c(9, 6, 6, 3)
############################

##############SAVING FILES##############
write.csv(barcodes.refalt, file = "out/barcodes-refalt.csv", quote = F, row.names = F)
write.csv(barcodes.refalt.filtered, file = "out/barcodes-refalt-filtered.csv", quote = F, row.names = F)
write.csv(barcodes, file = "out/barcodes.csv", quote = F, row.names = F)
write.csv(barcodes.original, file = "out/barcodes-original.csv", quote = F, row.names = F)
snpavailable <- colnames(barcodes[startPosindex:endPosIndex])
chnames <- sapply(snpavailable, FUN = function(x) {unlist(strsplit(x, split = "['.']"))[1]})
chpos <- sapply(snpavailable, FUN = function(x) {unlist(strsplit(x, split = "['.']"))[2]})
snpavailable <- cbind(chnames, chpos)
write.table(snpavailable, file = "out/SNP-available.tsv", sep="\t", row.names = F, quote = F, col.names = F)

####################################
  
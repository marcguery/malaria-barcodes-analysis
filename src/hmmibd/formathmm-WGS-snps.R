################WGS snps#################
bigbarcodes <- read.csv("../read/out/WGS-19262snps_barcodes.csv", 
                        header = T)

startPosindex <- 2
endPosIndex <- ncol(bigbarcodes)-2
bigbarcodes$Barcode <- apply(bigbarcodes[,c(startPosindex:endPosIndex)], MARGIN = 1,
                          FUN = function(x) {paste(x, sep="", collapse = "")})
################################

################WGS snps MAPPING#################
bigWGSmapping <- read.csv("../read/out/WGS-19262snps_mapping.csv")
colnames(bigWGSmapping) <- paste(colnames(bigWGSmapping), "bigWGS", sep=".")
################################

################MERGING MAPPING################
bigbarcodes <- merge(bigbarcodes, bigWGSmapping[,c(1,3)], 
                  by.x = "Internal.Sample.ID", by.y = "sample.bigWGS", all.x = T)

bigbarcodes.refalt <- bigbarcodes[which(bigbarcodes$Internal.Sample.ID%in%c("ALT.WGS", "REF.WGS")),]
bigbarcodes <- bigbarcodes[-which(bigbarcodes$Internal.Sample.ID%in%c("ALT.WGS", "REF.WGS")),]
################################

##############FILTERING LOW QUALITY BARCODES##############
goodsnps.bigbarcodes <- sapply(X = bigbarcodes$Barcode,
                   FUN = function(x) {nchar(gsub(pattern = "[X|-]", replacement = "", x = x))})
min(goodsnps.bigbarcodes) #This should be greater than 0
ctoffsnps.bigbarcodes <- 4001
#Distribution of number of good quality SNPs in barcodes
quantile(goodsnps.bigbarcodes)
#Percentage of barcodes having no good SNPs at all
ecdf(goodsnps.bigbarcodes)(ctoffsnps.bigbarcodes-1)*100
#number of barcodes having no good SNPs at all
length(which(goodsnps.bigbarcodes<=ctoffsnps.bigbarcodes))

bigbarcodes <- bigbarcodes[goodsnps.bigbarcodes >= ctoffsnps.bigbarcodes,]

####################################

#####FOR WGS SNPS#########

#FOR REF AND ALT
hmmibdbigrefalt <- t(bigbarcodes.refalt)

colnames(hmmibdbigrefalt) <- hmmibdbigrefalt[1,]
hmmibdbigrefalt <- data.frame(hmmibdbigrefalt[startPosindex:endPosIndex,])
hmmibdbigrefalt$chrom <- unlist(lapply(rownames(hmmibdbigrefalt),
                                    FUN = function(x) { unlist(strsplit(x, split = "[.]"))[1] }))
hmmibdbigrefalt$pos <- unlist(lapply(rownames(hmmibdbigrefalt),
                                  FUN = function(x) { unlist(strsplit(x, split = "[.]"))[2] }))

hmmibdbigrefalt$chrom <- sub(pattern = "Pf3D7_(0)?", replacement = "", hmmibdbigrefalt$chrom)
hmmibdbigrefalt$chrom <- sub(pattern = "_v3", replacement = "", hmmibdbigrefalt$chrom)

hmmibdbigrefalt$chrom <- as.integer(hmmibdbigrefalt$chrom)
hmmibdbigrefalt$pos <- as.integer(hmmibdbigrefalt$pos)

#These SNPs selected for their quality 
#will be compared with hmmIBD output of consensus barcodes
hmmibdwgssnps <- t(bigbarcodes)
colnames(hmmibdwgssnps) <- hmmibdwgssnps[1,]
hmmibdwgssnps <- data.frame(hmmibdwgssnps[startPosindex:endPosIndex,])
hmmibdwgssnps$chrom <- unlist(lapply(rownames(hmmibdwgssnps),
                                      FUN = function(x) { unlist(strsplit(x, split = "[.]"))[1] }))
hmmibdwgssnps$pos <- unlist(lapply(rownames(hmmibdwgssnps),
                                    FUN = function(x) { unlist(strsplit(x, split = "[.]"))[2] }))

hmmibdwgssnps$chrom <- sub(pattern = "Pf3D7_(0)?", replacement = "", hmmibdwgssnps$chrom)
hmmibdwgssnps$chrom <- sub(pattern = "_v3", replacement = "", hmmibdwgssnps$chrom)

hmmibdwgssnps$chrom <- as.integer(hmmibdwgssnps$chrom)
hmmibdwgssnps$pos <- as.integer(hmmibdwgssnps$pos)

hmmibdwgssnps <- cbind(hmmibdwgssnps[,c((ncol(hmmibdwgssnps)-1):ncol(hmmibdwgssnps))],
                        hmmibdwgssnps[,-c((ncol(hmmibdwgssnps)-1):ncol(hmmibdwgssnps))])

hmmibdwgssnps[,-c(1:2)] <- apply(X = hmmibdwgssnps[,-c(1:2)], MARGIN = 2,
                                  FUN = cv4ibd, ref=hmmibdbigrefalt$REF, alt=hmmibdbigrefalt$ALT)

hmmibdwgssnps <- hmmibdwgssnps[with(hmmibdwgssnps, order(chrom, pos)), ]

#Removal of one sample 3D7 related + 1 sample with bad format
hmmibdwgssnps <- hmmibdwgssnps[,c(1:2,2+which(!sub(".WGS", "", colnames(hmmibdwgssnps[,3:ncol(hmmibdwgssnps)]))
                                                   %in%sub(".WGS", "", samplestoremove$ID)))]
##############################

############WRITE DATA############
write.table(hmmibdwgssnps, file = "out/hmmIBD-WGS-snps.txt", sep="\t", quote = F, col.names = T, row.names = F)
########################

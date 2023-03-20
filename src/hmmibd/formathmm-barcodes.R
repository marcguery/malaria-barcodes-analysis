barcodes.refalt <- read.csv("../read/out/barcodes-refalt.csv")
barcodes <- read.csv("../read/out/barcodes-refalt-filtered.csv")
samplestoremove <- read.table("../read/rawdata/samples-to-remove.txt", sep ="\t", header = TRUE)

#startPosindex and endPosIndex delimit the locations of the columns storing
#the bases called by SNP genotyping or WGS
startPosindex <- 2
endPosIndex <- ncol(barcodes)-6
######################CONVERT FORMAT FOR HMMIBD###############
#All positions that are unknown or mixed are ignored
cv4ibd <- function(sample, ref, alt){
  sample <- as.vector(sample)
  sample[sample==ref] <- 1
  sample[sample==alt] <- 0
  sample[sample=="N"] <- -1
  sample[sample=="-"] <- -1
  sample[sample=="X"] <- -1
  return(sample)
}

#FOR REF AND ALT
hmmibdrefalt <- t(barcodes.refalt)

colnames(hmmibdrefalt) <- hmmibdrefalt[1,]
hmmibdrefalt <- data.frame(hmmibdrefalt[startPosindex:endPosIndex,])
hmmibdrefalt$chrom <- unlist(lapply(rownames(hmmibdrefalt),
                                    FUN = function(x) { unlist(strsplit(x, split = "[.]"))[1] }))
hmmibdrefalt$pos <- unlist(lapply(rownames(hmmibdrefalt),
                                  FUN = function(x) { unlist(strsplit(x, split = "[.]"))[2] }))

hmmibdrefalt$chrom <- sub(pattern = "Pf3D7_(0)?", replacement = "", hmmibdrefalt$chrom)
hmmibdrefalt$chrom <- sub(pattern = "_v3", replacement = "", hmmibdrefalt$chrom)

hmmibdrefalt$chrom <- as.integer(hmmibdrefalt$chrom)
hmmibdrefalt$pos <- as.integer(hmmibdrefalt$pos)

#FOR BARCODES
hmmibdbarcodes <- t(barcodes)
colnames(hmmibdbarcodes) <- hmmibdbarcodes[1,]
hmmibdbarcodes <- data.frame(hmmibdbarcodes[startPosindex:endPosIndex,])
hmmibdbarcodes$chrom <- unlist(lapply(rownames(hmmibdbarcodes),
                                      FUN = function(x) { unlist(strsplit(x, split = "[.]"))[1] }))
hmmibdbarcodes$pos <- unlist(lapply(rownames(hmmibdbarcodes),
                                    FUN = function(x) { unlist(strsplit(x, split = "[.]"))[2] }))

hmmibdbarcodes$chrom <- sub(pattern = "Pf3D7_(0)?", replacement = "", hmmibdbarcodes$chrom)
hmmibdbarcodes$chrom <- sub(pattern = "_v3", replacement = "", hmmibdbarcodes$chrom)

hmmibdbarcodes$chrom <- as.integer(hmmibdbarcodes$chrom)
hmmibdbarcodes$pos <- as.integer(hmmibdbarcodes$pos)

hmmibdbarcodes <- cbind(hmmibdbarcodes[,c((ncol(hmmibdbarcodes)-1):ncol(hmmibdbarcodes))],
                        hmmibdbarcodes[,-c((ncol(hmmibdbarcodes)-1):ncol(hmmibdbarcodes))])

hmmibdbarcodes[,-c(1:2)] <- apply(X = hmmibdbarcodes[,-c(1:2)], MARGIN = 2,
                                  FUN = cv4ibd, ref=hmmibdrefalt$REF, alt=hmmibdrefalt$ALT)

hmmibdbarcodes <- hmmibdbarcodes[with(hmmibdbarcodes, order(chrom, pos)), ]

hmmibdbarcodes.usedsamples <- hmmibdbarcodes[,c(1:2,2+which(!colnames(hmmibdbarcodes[,c(3:ncol(hmmibdbarcodes))])%in%samplestoremove$ID))]
########################

############WRITE DATA############
write.table(hmmibdbarcodes.usedsamples, file = "out/hmmIBD-barcodes-filtered.txt", sep="\t", quote = F, col.names = T, row.names = F)
write.table(hmmibdbarcodes, file = "out/hmmIBD-barcodes-refalt.txt", sep="\t", quote = F, col.names = T, row.names = F)
########################

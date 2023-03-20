######################READING######################
if (datatype == "WGS"){
  barcodes <- read.csv("../read/out/WGS-19262snps_barcodes.csv")
  mapping <- read.csv("../read/out/WGS-19262snps_mapping.csv")
  mapping <- mapping[!is.na(mapping$SampleID_corrected),]
  
  barcodes.ref <- barcodes[which(barcodes$Internal.Sample.ID%in%c("ALT.WGS", "REF.WGS")),]
  rownames(barcodes.ref) <- barcodes.ref$Internal.Sample.ID
  barcodes.ref <- barcodes.ref[,-c(1, (ncol(barcodes.ref)-2):ncol(barcodes.ref))]
  
  barcodes <- merge(barcodes, mapping[,c("sample", "SampleID_corrected")], 
                     by.x = "Internal.Sample.ID",
                     by.y = "sample")
  barcodes <- barcodes[!barcodes$Internal.Sample.ID%in%samplestoremove$ID,]
  rownames(barcodes) <- barcodes$SampleID_corrected
  barcodes <- barcodes[,-c(1, (ncol(barcodes)-2):ncol(barcodes))]
  #To get only samples that have FWS value (they have > 4k SNPs)
  barcodes <- barcodes[row.names(barcodes)%in%Fws$Sample,]
  
  colnames(barcodes) <- sub("\\.", ":", colnames(barcodes))
  barcodes <- data.frame(t(barcodes))
  
  newcols <- sapply(colnames(snps.alt), 
                    function(x){
                      mapping$SampleID_corrected[sub("\\.WGS", "", mapping$sample) == x]})
  if (length(unique(newcols)) == ncol(snps.alt)){
    colnames(snps.alt) <- newcols
  }
  
  newcols <- sapply(colnames(snps.ref), 
                    function(x){
                      mapping$SampleID_corrected[sub("\\.WGS", "", mapping$sample) == x]})
  if (length(unique(newcols)) == ncol(snps.ref)){
    colnames(snps.ref) <- newcols
  }
  
  newcols <- sapply(colnames(snps.sum), 
                    function(x){
                      mapping$SampleID_corrected[sub("\\.WGS", "", mapping$sample) == x]})
  if (length(unique(newcols)) == ncol(snps.sum)){
    colnames(snps.sum) <- newcols
  }
  
  
  barcodes <- barcodes[row.names(snps.alt),colnames(snps.alt)]
  
}else if (datatype == "SNP"){
  barcodes <- read.csv("../read/out/barcodes-refalt-filtered.csv")
  
  barcodes.ref <- barcodes[which(barcodes$Internal.Sample.ID%in%c("ALT.WGS", "REF.WGS")),]
  rownames(barcodes.ref) <- barcodes.ref$Internal.Sample.ID
  barcodes.ref <- barcodes.ref[,-c(1, (ncol(barcodes.ref)-5):ncol(barcodes.ref))]
  
  barcodes <- barcodes[!barcodes$Internal.Sample.ID%in%samplestoremove$ID,]
  barcodes <- barcodes[!barcodes$Internal.Sample.I%in%c("ALT.WGS", "REF.WGS"),]
  rownames(barcodes) <- barcodes$ID
  barcodes <- barcodes[,-c(1, (ncol(barcodes)-5):ncol(barcodes))]
  
  colnames(barcodes) <- sub("\\.", ":", colnames(barcodes))
  barcodes <- data.frame(t(barcodes))
}else{
    print("Unrecognized datatype")
  }
############################################

######################UTILS######################
source("rawdata/THEREALMcCOIL/categorical_method/McCOIL_categorical.R")

############################################

################FILTER SNPS##########################
basetozyg <- function(barcode, alt, ref){
  barcode[which(barcode==alt)] <- 0
  barcode[which(barcode==ref)] <- 1
  barcode[which(barcode=="N")] <- 0.5
  barcode[which(barcode=="X")] <- -1
  return(barcode)
}
  
if (datatype == "WGS"){
  barcodes[c(1:nrow(barcodes)), c(1:ncol(barcodes))] <- -1
  #Format for mccoil
  # 0.5 = mixed, 0 = maj allele, 1 = min allele, -1 = NA
  majoraf <- pmax(snps.alt, snps.ref)/(snps.ref+snps.alt)
  barcodes[majoraf < 0.95] <- 0.5
  barcodes[snps.ref > snps.alt & snps.sum > lowdepth & majoraf >= 0.95] <- 0
  barcodes[snps.alt > snps.ref & snps.sum > lowdepth & majoraf >= 0.95] <- 1
  
}else if (datatype == "SNP"){
  #Format for mccoil
  # 0.5 = mixed, 0 = maj allele, 1 = min allele, -1 = NA
  barcodes <- data.frame(apply(barcodes,
                      2, basetozyg,
                      barcodes.ref[which(rownames(barcodes.ref)=="ALT.WGS"),],
                      barcodes.ref[which(rownames(barcodes.ref)=="REF.WGS"),]))
  barcodes <- data.frame(apply(barcodes, c(1,2), as.numeric))
}else{
  print("Unrecognized datatype")
}


#Remove SNPs if present in less than 80% of samples
leastunk <- apply(barcodes, 1, function(x){
  length(which(x == -1))/length(x) < 0.2
  })
#Remove SNPs if Population MAF is > 0.25
# AF = reads (REF or ALT) + 1/2 reads (REF+ALT)
# mostvar <- apply(barcodes, 1, function(x){
#   abs((length(which(x == 0))+
#          0.5*length(which(x == 0.5)))
#       /(length(which(x == 0)) +
#           length(which(x == 1)) + 
#           0.5*length(which(x == 0.5))) - 0.5) < 0.25
# })
barcodes <- barcodes[leastunk,]

#Remove SNPs if in a window of 5kb (to avoid LD)
newbarcodedata <- barcodes
newbarcodedata = newbarcodedata[F,]
for (chr in unique(sub(":\\S+", "", row.names(barcodes)))){
  print(chr)
  barcodes.subset <- barcodes[grepl(chr, row.names(barcodes)),]
  pos <- as.numeric(sub("\\S+:", "", row.names(barcodes.subset)))
  DM = as.matrix(dist(pos))
  diag(DM) = 5000            ## ignore diagonal
  filtrows <- which(DM < 5000, arr.ind=TRUE)
  if(nrow(filtrows) > 0){
    print(head(filtrows))
    barcodes.subset <- barcodes.subset[-c(filtrows[filtrows[,1] < filtrows[,2],2]),]
  }
  newbarcodedata <- rbind(newbarcodedata, barcodes.subset)
}
barcodes <- newbarcodedata
rm(newbarcodedata)

barcodes <- t(barcodes)
#mcCOIL run
#Even with high number of runs (up to 30k) or burnin (up to 2k)
#some samples do not converge (sd too high)
#Results tend to not vary with initial COI
McCOIL_categorical(barcodes,
                   maxCOI=30, 
                   threshold_ind=0, threshold_site=0, 
                   totalrun=10000, burnin=1000, 
                   M0=15, e1=0.01, e2=0.01, 
                   path=paste0(getwd(), "/rawdata/THEREALMcCOIL/categorical_method"), 
                   output=paste0("../../../out/", datatype, "-mccoil_output.txt" ))
############################################
##################INTRO##################
#Generates the WGS barcodes
####################################
##################READ##################
snpsfile <- "rawdata/WGS-101snps"
outfile <- "WGS-barcodes"
snps <- read.table(paste(snpsfile, ".vcf", sep=""), 
                   sep="\t", header = T, skip = 52, comment.char = "")
snps.ad <- read.table(paste(snpsfile, "_allele-depth.tsv", sep=""), 
                      sep="\t",header = T)
####################################

##################SHAPE##################
snps <- snps[,c(1:5)]
colnames(snps)[1] <- "CHROM"

snps.merged <- merge(snps, snps.ad, by = c("CHROM", "POS"))
####################################

##################OPTIONS##################
lowdepth <- 6 #if depth is inferior strict to NUM, SNP is 'X'
if (!exists("minprop")){
  minprop <- 0.8
  } #if prop is inferior strict to NUM, SNP is 'N'
#minprop is set up from the pipeline calling this script
####################################
##################FUNCTIONS##################
parsedepth <- function(ad){
  intlist <- as.integer(unlist(strsplit(ad, ",")))
  res <- intlist[1]+intlist[2]
  return(c(intlist[1], intlist[2], res))
}

getgenotype <- function(ad){
  depth <- parsedepth(ad)
  refdepth <- depth[1]
  altdepth <- depth[2]
  totdepth <- depth[3]
  if (max(refdepth, altdepth) < lowdepth){
    return("X")
  }
  if (abs((refdepth/totdepth)-0.5) < (minprop-0.5)){
    return("N")
  }
  if (refdepth >= altdepth){
    return('ref')
  }else{
    return('alt')
  }
}

convertbarcode <- function(barcode, ref, alt){
  barcode[barcode=="ref"] <- ref[barcode=="ref"]
  barcode[barcode=="alt"] <- alt[barcode=="alt"]
  
  return(barcode)
}
####################################
##################USING AD TO GET GENOTYPE##################
gbarcodes <- snps.merged
gbarcodes[,-c(1:5)] <- apply(snps.merged[,-c(1:5)], MARGIN = c(1,2),
               FUN = getgenotype)
#Conversion
gbarcodes[,-c(1:5)] <- apply(gbarcodes[,-c(1:5)], MARGIN = 2,
                   FUN = convertbarcode, 
                   ref=as.vector(snps.merged$REF), alt=as.vector(snps.merged$ALT))
#################TRANSFORMING DATAFRAME#################
sample.gbarcodes <- t(gbarcodes)
colnames(sample.gbarcodes) <- paste(sample.gbarcodes[1,], 
                                    sub(pattern =" ", replacement = "", sample.gbarcodes[2,]), 
                                        sep = ".")
sample.gbarcodes <- as.data.frame(sample.gbarcodes[-c(1,2,3),])
sample.gbarcodes <- cbind(row.names(sample.gbarcodes), sample.gbarcodes)
colnames(sample.gbarcodes)[1] = "Internal.Sample.ID"

sample.gbarcodes$Barcode <- apply(X = sample.gbarcodes[,-1], MARGIN = 1,
                                  FUN =  function(x) {paste(x, sep="", collapse = "")})
sample.gbarcodes$Study <- "WGS"
########################################################

#########################MAPPING#########################
mapping <- read.csv("rawdata/WGS-snps_mapping.csv", header = T)[,c(1,3,4)]

#Line shifting correction
barcodemapping <- read.csv("out/SNP-barcodes_mapping.csv")
barcodemapping <- barcodemapping[c(351:576),]

sameIDs <- merge(mapping, barcodemapping,
            by.x = "SampleID_corrected", by.y = "Full_ID")
nmapid <- which(barcodemapping$Full_ID%in%sameIDs$SampleID_corrected)
gmapid <- which(mapping$SampleID_corrected%in%sameIDs$SampleID_corrected)
mapping[gmapid, c("SampleID_corrected")] <- barcodemapping$Full_ID[nmapid+1]

#Removing bad formatted IDs
mapping$SampleID_corrected <- sub("16PCD", "1611", mapping$SampleID_corrected)
mapping$SampleID_corrected[which(mapping$SampleID_corrected=="j0111201_1412")] <- "J0111201_1412"
mapping$SampleID_corrected[which(mapping$SampleID_corrected=="K0344022_0000")] <- NA
#Replacing one wrongly assigned barcode with another
mapping$SampleID_corrected[which(mapping$SampleID_corrected=="K0262922_1704")] <- "K0262940_1704.2"
mapping$SampleID_corrected[which(mapping$SampleID_corrected=="K0262940_1704")] <- "K0262922_1704"
mapping$SampleID_corrected[which(mapping$SampleID_corrected=="K0262940_1704.2")] <- "K0262940_1704"
mapping$SampleID_corrected[!grepl(pattern = "[A-Z]\\d{7}_\\d{4}$", mapping$SampleID_corrected)] <- NA
####################################################
sample.gbarcodes$Internal.Sample.ID <- paste(sample.gbarcodes$Internal.Sample.ID, ".WGS", sep = "")
mapping$sample <- paste(make.names(mapping$sample), ".WGS", sep = "")
############SAVING FILES############
write.csv(sample.gbarcodes, file = paste("out/", outfile, ".csv", sep=""),
          quote=F, row.names = F)

write.csv(mapping, file = paste("out/", outfile, "_mapping.csv", sep=""), quote = F, row.names = F)
####################################


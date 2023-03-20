#This script needs functions from WGS-barcodes.R
############WGS SNPS##############
bigsnpsfile <- "WGS-19262snps"
headerlength <- 52
bigsnps.data <- read_delim(paste0("rawdata/", bigsnpsfile, ".vcf.zip"),
                           delim="\t", comment = "", skip = headerlength, 
                           col_types = cols(.default = "c") )

bigsnps.ad <- read_delim(paste0("rawdata/", bigsnpsfile, "_allele-depth.tsv.zip"),
                         delim="\t", col_types = cols(.default = "c"))

##################SHAPE##################
bigsnps.data <- bigsnps.data[,c(1:5)]
colnames(bigsnps.data)[1] <- "CHROM"

bigsnps.merged <- merge(bigsnps.data, bigsnps.ad, by = c("CHROM", "POS"), sort = F)
####################################

##################USING AD TO GET GENOTYPE##################
bigsnps <- bigsnps.merged
bigsnps[,-c(1:5)] <- apply(bigsnps.merged[,-c(1:5)], MARGIN = c(1,2),
                           FUN = getgenotype)
#Conversion
bigsnps[,-c(1:5)] <- apply(bigsnps[,-c(1:5)], MARGIN = 2,
                           FUN = convertbarcode, 
                           ref=as.vector(bigsnps.merged$REF), alt=as.vector(bigsnps.merged$ALT))
####################################

#################TRANSFORMING DATAFRAME#################
sample.bigsnps <- t(bigsnps)
colnames(sample.bigsnps) <- paste(sample.bigsnps[1,], 
                                  sub(pattern =" ", replacement = "", sample.bigsnps[2,]), 
                                  sep = ".")
sample.bigsnps <- as.data.frame(sample.bigsnps[-c(1,2,3),])
sample.bigsnps <- cbind(row.names(sample.bigsnps), sample.bigsnps)
colnames(sample.bigsnps)[1] = "Internal.Sample.ID"

sample.bigsnps$Barcode <- apply(X = sample.bigsnps[,-1], MARGIN = 1,
                                FUN =  function(x) {paste(x, sep="", collapse = "")})
sample.bigsnps$Study <- "BigWGS"
########################################################

###############REMOVING APICOPLAST CHROMOSOME###############
#These samples are SNPs from apicolplast
grep("API", colnames(sample.bigsnps))
sample.bigsnps <- sample.bigsnps[,-c(grep("API", colnames(sample.bigsnps)))]
##############################

sample.bigsnps$Internal.Sample.ID <- make.names(sample.bigsnps$Internal.Sample.ID)
sample.bigsnps$Internal.Sample.ID <- paste(sample.bigsnps$Internal.Sample.ID, ".WGS", sep = "")

############SAVING FILES############
write.csv(sample.bigsnps, file = paste("out/", bigsnpsfile, "_barcodes.csv", sep=""),
          quote=F, row.names = F)
write.csv(mapping, file = paste("out/", bigsnpsfile, "_mapping.csv", sep=""),
          quote = F, row.names = F)
####################################
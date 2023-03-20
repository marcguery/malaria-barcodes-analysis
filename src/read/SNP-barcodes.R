###############BARCODES###############
#Loading raw data files
options(stringsAsFactors = F)
barcodefiles <- dir(path = "rawdata", pattern = "SNP-barcodes[1-3]{1}.csv")
barcode.1 <- read.csv(paste0("rawdata/",barcodefiles[1]))
barcode.2 <- read.csv(paste0("rawdata/",barcodefiles[2]))
barcode.3 <- read.csv(paste0("rawdata/",barcodefiles[3]))

#Removing the first column containing the study name
barcode.1 <- barcode.1[,-1]
barcode.1$Study <- "SNP"
barcode.2$Study <- "SNP"
barcode.3$Study <- "SNP"
#Updating the column name to match the other studies
names(barcode.2)[names(barcode.2)=="ID"] <- "Internal.Sample.ID"

barcodes <- rbind(barcode.1, barcode.2)
barcodes <- rbind(barcodes, barcode.3)
barcodes[is.na(barcodes)] <- "X"
#Building a continuous string to store barcodes
barcodes$Barcode <- apply(barcodes[,-c(1,2,ncol(barcodes))], 1, function(x){paste(x, collapse = "")})

#Checking if a character in a barcode is unexpected
unexpectedchar <- apply(X = barcodes[,-which(colnames(barcodes)%in%c("Internal.Sample.ID", "Barcode", "Study"))], 
                        FUN = function (x) { x!="-" & x!="X" & x!="N" & x!="A" & x!="T" & x!="G" & x!="C" }, 
                        MARGIN = 2)
boolunexpectedchar <- apply(X = unexpectedchar, 
                            FUN = function(x) {any(x)}, 
                            MARGIN = 1)
which(boolunexpectedchar) #it should be empty
##############################

###############ID MAPPING###############
#Loading raw data files
mappingfiles <- dir(path = "rawdata", pattern = "SNP-barcodes_mapping[1-3]{1}.csv")
mapping.1 <- read.csv(paste0("rawdata/",mappingfiles[1]))[,c(2,3,5)]
mapping.2 <- read.csv(paste0("rawdata/",mappingfiles[2]))[,c(2,3,4)]
mapping.3 <- read.csv(paste0("rawdata/",mappingfiles[3]))[,c(2,3,5)]

colnames(mapping.2) <- colnames(mapping.1)

######LINE shifting
#From SPT14589 to SPT14814, there could have had a line skip...
firstbad <- which(mapping.1$Internal.Sample.ID=="SPT14589")#Line 352 in the excel file
lastbad <- which(mapping.1$Internal.Sample.ID=="SPT14814")#Line 577 in the excel file
shifttable <- mapping.1[firstbad:lastbad,]
mapping.1$Full_ID[firstbad:lastbad] <- c(mapping.1$Full_ID[(firstbad+1):lastbad], NA)
shifttable$Full_ID_after_shift <- mapping.1$Full_ID[firstbad:lastbad]
colnames(shifttable)[3] <- "Full_ID_before_shift"
write.csv(shifttable, file="out/SNP-barcodes_old-vs-new-shifted-IDs.csv", quote = F, row.names = F)
write.csv(mapping.1, file = "out/SNP-barcodes_mapping1-shifted.csv", quote=F, row.names = F)
######

mappings <- rbind(mapping.1, mapping.3)
mappings <- rbind(mappings, mapping.2)

#Removing bad formatted IDs
mappings$Full_ID[grepl("3D7", mappings$Full_ID)] <- NA
mappings$Full_ID[grepl("MIX", mappings$Full_ID)] <- NA
mappings$Full_ID[grepl("/", mappings$Full_ID)] <- NA
mappings$Full_ID <- sub("16PCD", "1611", mappings$Full_ID)
mappings$Full_ID <- sub("1607F", "1607", mappings$Full_ID)
mappings$Full_ID <- sub("1610F", "1610", mappings$Full_ID)
mappings$Full_ID <- sub("a_", "_", mappings$Full_ID)
mappings$Full_ID <- sub("b_", "_", mappings$Full_ID)
mappings$Full_ID[which(mappings$Full_ID=="j0111201_1412")] <- "J0111201_1412"
mappings$Full_ID[which(mappings$Full_ID=="KO090903_1607")] <- "K0090903_1607"
mappings$Full_ID[!grepl(pattern = "[A-Z]\\d{7}_\\d{4}$", mappings$Full_ID)] <- NA
#Replacing one wrongly assigned barcode with another
mappings$Full_ID[which(mappings$Full_ID=="K0262922_1704")] <- "K0262940_1704.2"
mappings$Full_ID[which(mappings$Full_ID=="K0262940_1704")] <- "K0262922_1704"
mappings$Full_ID[which(mappings$Full_ID=="K0262940_1704.2")] <- "K0262940_1704"
#Number of blood samples
nrow(mappings[grepl("[A-z]\\d{7}_1\\d{1,3}", mappings$Full_ID),])
#Number of participants
length(unique(sub("_\\S+", "", mappings$Full_ID[grepl("[A-z]\\d{7}_1\\d{1,3}", mappings$Full_ID)])))

#Showing duplicated IDs
dupids <- which(duplicated(mappings$Full_ID, fromLast = F))
dupids
#30 146 282 297 310 318 358 499 564

###############################################
##################SAVING FILES##############
all(barcodes$Internal.Sample.ID%in%mappings$Internal.Sample.ID)
mappings[which(!mappings$Internal.Sample.ID%in%barcodes$Internal.Sample.ID),]

write.csv(barcodes, "out/SNP-barcodes.csv",
          quote = F, row.names = F) #Barcodes without any filter

write.csv(x = mappings, file = "out/SNP-barcodes_mapping.csv", quote = F, row.names = F)

###############################################

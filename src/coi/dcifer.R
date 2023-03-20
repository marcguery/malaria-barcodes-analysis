###################READING###################
locusdata.snp <- read.table("../hmmibd/out/hmmIBD-WGS-snps.txt", h = T)

locusdata.locus <- data.frame(t(data.frame(apply(locusdata.snp[,-c(1,2)], 1, 
                         FUN = function(x){
                           freqs = c(length(which(x == 0)),
                                     length(which(x == 1)),
                                     length(which(x == -1)))
                           return(c("ref" = freqs[1], "alt" = freqs[2], "unk" = freqs[3]))
                         }))))

leastunk <- locusdata.locus$unk/(locusdata.locus$ref + locusdata.locus$alt  + locusdata.locus$unk) < 0.5
mostvar <- abs(locusdata.locus$ref/(locusdata.locus$ref + locusdata.locus$alt) - 0.5) < 0.2
locusdata.locus <- locusdata.locus[leastunk & mostvar,]
locusdata.snp <- locusdata.snp[which(leastunk & mostvar),]

newlocusdata.snp <- locusdata.snp
newlocusdata.snp <- newlocusdata.snp[F,]
for (chr in unique(locusdata.snp$chrom)){
  print(chr)
  locusdata.snp.subset <- locusdata.snp[locusdata.snp$chrom == chr,]
  DM = as.matrix(dist(locusdata.snp.subset$pos))
  diag(DM) = 1000            ## ignore diagonal
  filtrows <- which(DM < 1000, arr.ind=TRUE)
  locusdata.snp.subset <- locusdata.snp.subset[-c(filtrows[filtrows[,1] < filtrows[,2],2]),]
  newlocusdata.snp <- rbind(newlocusdata.snp, locusdata.snp.subset)
}
locusdata.snp <- newlocusdata.snp
rm(newlocusdata.snp)


locusdata.snp.melt <- melt(locusdata.snp, id.vars = c("chrom", "pos"), 
                           variable.name = "sample",
                           value.name = "allele")
locusdata.snp.melt$locus <- paste(locusdata.snp.melt$chrom, locusdata.snp.melt$pos, sep = ".")
locusdata.snp.melt <- locusdata.snp.melt[,c("sample", "locus", "allele")]
locusdata.snp.melt <- locusdata.snp.melt[locusdata.snp.melt$allele != -1,]
######################################
###################LOADING###################
dsmp <- formatDat(locusdata.snp.melt, svar = "sample", lvar = "locus", avar = "allele")
str(dsmp, list.len = 2)
######################################
###################COI AND AF###################
lrank <- 2
coi   <- getCOI(dsmp, lrank = lrank)
afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 
str(afreq, list.len = 2)
######################################
###################RELATEDNESS###################
dres0 <- ibdDat(dsmp, coi, afreq, pval = TRUE, confint = TRUE, rnull = 0, 
                alpha = 0.05, nr = 1e3)
View(dres0)
######################################


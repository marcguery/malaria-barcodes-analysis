#!/usr/bin/env Rscript
#27/09/2021

#This script will scrutinize genetic diversity of malaria parasite Plasmodium falciparum
#in 4 nearby villages of the eastern part of The Gambia between 2014 and 2017

#Steps of the analysis pipeline:

#1. Read and homogenize data input files of epidemiological, barcode and WGS data

#2. Format barcode files to be read by hmmIBD and calculate IBDs

#3. Build networks of samples based on IBD on barcodes

#4. Retrieve strain genetic history with WGS 19k SNPs

#5. Estimate COI using barcode consensus 101 SNPs and WGS 19k SNPs

#6. Correlate IBD between barcode consensus 101 SNPs and WGS 19k SNPs

library("reshape2")
library("lubridate")
library("cowplot")
library("scales")
library("dplyr")
library("ggplot2")
library("ggnewscale")
library("ggforce")
library("readr")
library("stringr")
library("reticulate")
library("igraph")
library("gganimate")
library("ggtext")
library("IRanges")
library("RColorBrewer")
library("ggpattern")
#need also package "sets" which should not be loaded here
# as there are namespace conflicts

#need also package "GenomicRanges"
options(stringsAsFactors = F, tz="GMT")
setwd("./src/read")
###########################1. READ###########################
setwd("../read")
source("SNP-barcodes.R") #Load genotyped barcodes
source("WGS-barcode-mixed-cutoff.R") #Estimate WGS barcodes N cutoff
source("WGS-barcodes.R") #Load WGS snp for barcodes
source("WGS-SNPs.R") #Load WGS whole good quality SNPs
source("merge-SNP-WGS.R") #Merge barcodes from genotyping and WGS
source("figures.R")
rm(list=ls())
######################################################
###########################2. HMMIBD###########################
setwd("../hmmibd")
source("formathmm-barcodes.R") #Format barcode data for hmmIBD
source("formathmm-WGS-snps.R") #Format WGS SNP data for hmmIBD
rm(list=ls())
#system("./run-hmmibd-local.sh", wait=FALSE)

######################################################
###########################3. NETWORKS###########################
setwd("../network")
source("barcodes-network.R") #Build network
source("diversity.R")  #Get the distribution of IBD values
source("network-filter.R") #Filter network
source("network-doi-strains.R") #Estimate individual and population 
                                #duration of infection
source("remove-continf.R")
removecontinf <- FALSE
source("group-network.R") #Group samples by location and date
source("figures.R")
rm(list=ls())
######################################################
###########################4. TRIADES###########################
setwd("../triades")
parentalmin <- 0.35
parentalmax <- 0.65
identical <- 0.9
different <- 0.2
source("WGS-snps-networks.R") #Create network from WGS SNPs
source("print-IBD.R") #Find triades from WGS SNPs
source("find-triades.R") #Find triades from WGS SNPs
source("show-multiads.R") #Show IBD fragments shared by 2+ isolates
source("get-IBD-matrtix.R") #Summary of IBD shared fragments
rm(list = ls())
######################################################

###########################5. COI###########################
setwd("../coi")
source("fws.R")
datatype <- "SNP"
source("mccoil.R")
datatype <- "WGS"
source("mccoil.R")
source("figures.R")
rm(list = ls())
######################################################

###########################6. CORRELATION###########################
setwd("../ibdcorrelation")
source("WGS-barcodes-correlation.R") #Compare IBDs between genotyped and WGS samples
rm(list = ls())
######################################################

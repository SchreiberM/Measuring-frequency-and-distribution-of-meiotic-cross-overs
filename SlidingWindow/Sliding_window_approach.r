# Author: Miriam Schreiber
# Script for identification of recombination breakpoints
# Needs a .tsv input

# Provide necessary libraries
library(ggplot2)
library(zoo)
library(dplyr)
library(BBmisc)
library(runner)
library(spatstat)


# This script is adjusted to work with barley sequencing data mapped to Morex V2 in split chromosome format
# The first part of the script merges the split chromosome information into one for each of the seven chromosomes

# Read in file
# Each file will consist of lines with the following information
# chromosome    position    heterozygosity (0.5 for heterozygous or 0 for homozygous)   colour (grey for homozygous, blue for heterozygous)
# for example: chr1H_part1	1617271	0.5	blue

# snps <- read.table("YYD5.filtered_2.tsv",sep="\t",header=F,stringsAsFactors=FALSE)
snps <- read.table(<file>,sep="\t",header=F,stringsAsFactors=FALSE)

# Provide sample name
# name <- "Bowman_EMS_M3_1"
name <- <name>


#########################################

#set the column namname <- "Bowman_EMS_M3_3"
colnames(snps)<-c("chr","start","GENO","COLOUR")

snps$chr <- as.character(snps$chr)
snps$GENO <- as.double(snps$GENO)

# chromosomes have been split into two parts, combine into one (here it is using Morex V2)
snps$chr[snps$chr == "chr1H_part1"] <- "chr1H"
snps$start[snps$chr == "chr1H_part2"] <- snps$start[snps$chr == "chr1H_part2"] + 312404992
snps$chr[snps$chr == "chr1H_part2"] <- "chr1H"

snps$chr[snps$chr == "chr2H_part1"] <- "chr2H"
snps$start[snps$chr == "chr2H_part2"] <- snps$start[snps$chr == "chr2H_part2"] + 341635707
snps$chr[snps$chr == "chr2H_part2"] <- "chr2H"

snps$chr[snps$chr == "chr3H_part1"] <- "chr3H"
snps$start[snps$chr == "chr3H_part2"] <- snps$start[snps$chr == "chr3H_part2"] + 347769808
snps$chr[snps$chr == "chr3H_part2"] <- "chr3H"

snps$chr[snps$chr == "chr4H_part1"] <- "chr4H"
snps$start[snps$chr == "chr4H_part2"] <- snps$start[snps$chr == "chr4H_part2"] + 326411310
snps$chr[snps$chr == "chr4H_part2"] <- "chr4H"

snps$chr[snps$chr == "chr5H_part1"] <- "chr5H"
snps$start[snps$chr == "chr5H_part2"] <- snps$start[snps$chr == "chr5H_part2"] + 317707318
snps$chr[snps$chr == "chr5H_part2"] <- "chr5H"

snps$chr[snps$chr == "chr6H_part1"] <- "chr6H"
snps$start[snps$chr == "chr6H_part2"] <- snps$start[snps$chr == "chr6H_part2"] + 314597075
snps$chr[snps$chr == "chr6H_part2"] <- "chr6H"

snps$chr[snps$chr == "chr7H_part1"] <- "chr7H"
snps$start[snps$chr == "chr7H_part2"] <- snps$start[snps$chr == "chr7H_part2"] + 336913125
snps$chr[snps$chr == "chr7H_part2"] <- "chr7H"


#order the chromosomes
goodChrOrder <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "chrUn")

#apply this as a factor
snps$chr <- factor(snps$chr,levels=goodChrOrder)

# write("subset the data to exclude chrUn", stdout())
snps_subset <- snps[grep ("chrUn", snps$chr, invert=TRUE), ]
snps_subset2 <- na.omit(snps_subset)


#########################################

# The second part of the script is the sliding window approach. The SNP dataset is split into chunks (chunkSizes)
# Based on the number of SNPs within a chunk the bp windowsize is determined (BPSize)
# A weighted median is calculated across a sliding window using the runner function. 

newResults <- data.frame()

#sliding window appraoch to determine cross over positions
for(a in 1:length(unique(snps_subset$chr))){
    
    chromosome <- snps_subset2[snps_subset2$chr==unique(snps_subset$chr)[a],]
    chromosome_temp <- snps_subset2[snps_subset2$chr==unique(snps_subset$chr)[a],]

    # define chunk sizes for adjusted window size in Mb
    chunkSizes <- chunk(chromosome$start,props=c(8,12,18,24,18,12,8))
    chromosome$BPSize <- NA
    chromosome$BPSize[chromosome$start%in%chunkSizes[[1]]] <- length(chunkSizes[[1]])/0.00001
    chromosome$BPSize[chromosome$start%in%chunkSizes[[2]]] <- length(chunkSizes[[2]])/0.00001
    chromosome$BPSize[chromosome$start%in%chunkSizes[[3]]] <- length(chunkSizes[[3]])/0.00001
    chromosome$BPSize[chromosome$start%in%chunkSizes[[4]]] <- length(chunkSizes[[4]])/0.00001
    chromosome$BPSize[chromosome$start%in%chunkSizes[[5]]] <- length(chunkSizes[[5]])/0.00001
    chromosome$BPSize[chromosome$start%in%chunkSizes[[6]]] <- length(chunkSizes[[6]])/0.00001
    chromosome$BPSize[chromosome$start%in%chunkSizes[[7]]] <- length(chunkSizes[[7]])/0.00001
    
    # define weight for weighted median in function below
    chromosome$weight[chromosome$GENO==0] <- 1
    chromosome$weight[chromosome$GENO==0.5] <- 1.5

    adjustedGenoCall <- runner::runner(x=chromosome[,c(3,6)], k=chromosome$BPSize, lag=1, idx=chromosome$start, f=function(x) {weighted.median(x[,1], x[,2])})

    chromosome_temp$GENO <- adjustedGenoCall

    newResults <- rbind(newResults, chromosome_temp)
}


# output results are adjusted based on the new values.
newResults$GENO[newResults$GENO==0.25] <- 0.5
newResults$COLOUR[newResults$GENO==0] <- "grey"
newResults$COLOUR[newResults$GENO==0.5] <- "blue"

# For a check how many positions were changed by the sliding window approach a quick check can be done
# FALSE are the number of positions which are different to the raw data 
# table(snps_subset2$COLOUR == newResults$COLOUR)

#########################################

# Plot raw zygosity call of SNPs along the 7 barley chromosomes
piPlot_raw <- ggplot(snps_subset2, aes(start, GENO)) + 
geom_point(size = 0.4, colour=snps_subset2$COLOUR,aes(x=start, y=GENO)) + facet_wrap(~ chr,ncol=1) +
ylim(-0.5,1.5) +
ggtitle(paste(name, "raw variant call", sep=" ")) + 
scale_y_continuous(limits=c(-0.5,1),labels=c("", "Hom", "Het", "")) +
xlab("Physical position in bp") + ylab("") +
theme_classic(base_size = 8) 


tiff(paste(name, "rawVariant.tiff", sep="_"), units="cm", res=300, height=15, width=12, compression="lzw")
print(piPlot_raw)
dev.off()

# Plot the adjusted zygosity call of SNPs along the 7 barley chromosomes
piPlot<- ggplot(newResults, aes(start, GENO)) + 
geom_point(size = 0.4, colour=newResults$COLOUR,aes(x=start, y=GENO)) + facet_wrap(~ chr,ncol=1) +
ylim(-0.5,1.5) +
ggtitle(paste(name, "sliding_bp", sep=" ")) + 
xlab("Physical position in bp") + ylab("") +
scale_y_continuous(limits=c(-0.5,1),labels=c("", "Hom", "Het", "")) +
theme_classic(base_size = 8) 

tiff(paste(name, "sliding_bp.tiff", sep="_"), units="cm", res=300, height=15, width=12, compression="lzw")
print(piPlot)
dev.off()



#########################################
# Determine cross over positions
# The third part of the script uses the cleaned up results output to determine cross over positions

finalOut <- na.omit(newResults)

# If there is a dobule switch within 4 SNPs then merge togehter and do not count as cross over switch.

crossOverVector <- data.frame()

for(b in 1:7){
    
    chromosome <- finalOut[finalOut$chr==unique(finalOut$chr)[b],]

    switchPoint <- which(chromosome$COLOUR != dplyr::lag(chromosome$COLOUR))

    if(length(switchPoint)==0){
        temp <- cbind(as.character(chromosome$chr[1]),chromosome$start[1],chromosome$start[dim(chromosome)[1]],chromosome$GENO[dim(chromosome)[1]],chromosome$COLOUR[dim(chromosome)[1]])
        crossOverVector <- rbind(crossOverVector, temp)

    }else if(length(switchPoint)==1){
        temp <- cbind(as.character(chromosome$chr[1]),chromosome$start[1],chromosome$start[switchPoint-1],chromosome$GENO[switchPoint-1],chromosome$COLOUR[switchPoint-1])
        crossOverVector <- rbind(crossOverVector, temp)
        temp <- cbind(as.character(chromosome$chr[1]),chromosome$start[switchPoint],chromosome$start[dim(chromosome)[1]],chromosome$GENO[dim(chromosome)[1]],chromosome$COLOUR[dim(chromosome)[1]])
        crossOverVector <- rbind(crossOverVector, temp)

    }else{
        switchPointout <- switchPoint
        for(j in 2:length(switchPoint)){
            if((switchPoint[j]-switchPoint[j-1]) <= 4){
                if(is.na(switchPointout[j-1])){
                    switchPointout[j] <- switchPoint[j]
                }else if(length(switchPoint)==2){
                    switchPointout[j] <- NA
                }else{
                    switchPointout[j] <- NA
                    switchPointout[j-1] <- NA
                }
            }else if(switchPoint[j-1]==2){
                switchPointout[j-1] <- NA
            }
        }

        if((dim(chromosome)[1]-switchPoint[length(switchPoint)])<=3){
            switchPointout[length(switchPoint)] <- NA
        }

        switchPointCleaned <- na.omit(switchPointout)
        
        for(k in 1:(length(switchPointCleaned)+1)){
            if(k==1){
                temp <- cbind(as.character(chromosome$chr[1]),chromosome$start[1],chromosome$start[switchPointCleaned[k]-1],chromosome$GENO[switchPointCleaned[k]-1],chromosome$COLOUR[switchPointCleaned[k]-1])
                crossOverVector <- rbind(crossOverVector, temp)
            }else if(k==(length(switchPointCleaned)+1)){
                temp <- cbind(as.character(chromosome$chr[1]),chromosome$start[switchPointCleaned[k-1]],chromosome$start[dim(chromosome)[1]],chromosome$GENO[dim(chromosome)[1]],chromosome$COLOUR[dim(chromosome)[1]])
                crossOverVector <- rbind(crossOverVector, temp)
            }else{
                temp <- cbind(as.character(chromosome$chr[1]),chromosome$start[switchPointCleaned[k-1]],chromosome$start[switchPointCleaned[k]],chromosome$GENO[switchPointCleaned[k-1]],chromosome$COLOUR[switchPointCleaned[k-1]])
                crossOverVector <- rbind(crossOverVector, temp)
            }
        }
    }

}
colnames(crossOverVector) <- c("Chromosome", "Start", "End", "Value", "Colour")

# Count the number of cross over positions and prepare an output text file with the name and number of cross overs
RecNumber <- dim(crossOverVector)[1]-7

sink(paste(name, ".txt", sep=""))
print(name)

print(paste("The number of crossovers is:", RecNumber, sep=" "))

print(crossOverVector[,1:3])
sink()


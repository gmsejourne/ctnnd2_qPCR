library(xlsx)
library(stringr)
library(tidyverse)
library(ggplot2)
getwd()
setwd('/Users/gms50/Documents/Eroglu\ Lab/Data/qPCR/qPCR_220403')

##################
# Pre-processing #
##################
d <- read.xlsx("GMS_Lyst-geno-opt3_2022-04-03_162629.xls", 1, header=F)

dd <- d[44:(nrow(d)-4),] #-4 rows at the bottom that are info about the run
print(dd[1,1]) #check first row
print(dd[nrow(dd),1]) #check last row

x <- d[43, 1:ncol(d)]
colnames(dd) <- c()
  for (i in 1:ncol(dd)){
    colnames(dd)[i] <- as.character(x[[i]])
  }

print(colnames(dd[15])) #check column number of CT data
colnames(dd)[15] <- 'CT'
colnames(dd)[16] <- 'Ct Mean'
colnames(dd)[17] <- 'Ct SD'

print(dd$CT)
sapply(dd$CT, class) # check class that your data has defaulted to

for (i in 1:nrow(dd)){
  if (dd$CT[i] == 'Undetermined'){
    dd$CT[i] <- 10e10
  }
}
print(dd$CT)

for (i in 1:nrow(dd)){
  if (dd$`Ct Mean`[i] == ''){
    dd$`Ct Mean`[i] <- 10e10
  }
}
print(dd$`Ct Mean`)
  
dd$CT <- as.numeric(dd$CT) # convert data to numeric class
dd$`Ct Mean` <- as.numeric(dd$`Ct Mean`)
dd$`Ct SD` <- as.numeric(dd$`Ct SD`)
print(dd$CT) 
sapply(dd$CT, class) # check success of conversion

print(dd$`Sample Name`)
sapply(dd$`Sample Name`, class) # if Sample name is not already "character" class, convert using code below
#dd$`Sample Name` <- as.character(levels(dd$`Sample Name`))
#dd$`Target Name` <- as.character(levels(dd$`Target Name`))[dd$`Target Name`]

dd <- dd[!dd$Task == 'NTC', ] # omit NTC

# omit if over SD threshold
dd$omit <- c()
for (i in 1:nrow(dd)){
  dd$omit[i] <- ifelse(dd$`Ct SD`[i] > 0.8, 'yes', 'no')
}
print(dd$omit)
# dd <- dd[!dd$omit == 'yes', ] # omit high SD
print(dd$omit)

Dat <- aggregate(dd, by=list(dd$`Ct Mean`, dd$`Sample Name`, dd$`Target Name`), mean)
Dat <- Dat[, 1:3]

names(Dat) <- c('CtMean', 'Sample', 'Target')
print(Dat)

#############################################################
##                         Analysis                        ##
#############################################################

Dat$TargGroup <- rep(1:length(unique(Dat$Target)), each=length(unique(Dat$Sample)))
Dat$SampleGroup <- rep(1:length(unique(Dat$Sample)), times=length(unique(Dat$Target)))

# normalize to reference gene
k<-c()
for (i in 1:length(unique(Dat$TargGroup))){
  for (j in 1:length(unique(Dat$SampleGroup))){
    tmp <- Dat[Dat$SampleGroup == j, ]
    y <- tmp[tmp$TargGroup == i,1] - tmp[tmp$Target == 'ACTB', 1] #enter housekeeping gene
    k <- c(k, y)
  }
}
Dat$dCt <- k
rm(k,y, tmp)

# normalize to WT
k<-c()
for (i in 1:length(unique(Dat$TargGroup))){
  for (j in 1:length(unique(Dat$SampleGroup))){
    tmp <- Dat[Dat$TargGroup == i, ]
    y <- tmp[tmp$SampleGroup == j, 'dCt'] - tmp[tmp$SampleGroup == 1, 'dCt'] #enter sample to normalize to
    k <- c(k, y)
  }
}
Dat$ddCt <- k
Dat$RQ <- 2^(-Dat$ddCt)

write.csv(Dat, 'Lyst_genotyping_GMS_2022-04-03-analyzed.csv', row.names=F)

SampleGroup <- c(1:2)
Genotype <- c('WT', 'KO')
ad <- data.frame(SampleGroup, Genotype)
print(ad)

Dat <- merge(Dat, ad, by='SampleGroup')
print(Dat)

TargGroup <- c(1:7)
Exon <- c('B-actin', 'Exon 45', 'Exon 4-5', 'Exon 4-5', 'Exon 3-4', 'Exon 3-4', 'Exon 5')
bd <- data.frame(TargGroup, Exon)
print(bd)
Dat <- merge(Dat, bd, by='TargGroup')
print(Dat)

#p1 <- ggplot(Dat, aes(x=Genotype,y=RQ)) +
p1 <- ggplot(Dat[Dat$TargGroup > 1, ], aes(fill=Genotype,y=RQ,x=Exon)) + 
  geom_bar(stat='identity', position="dodge") + 
  #ylim(0, 1.2) + geom_jitter() + 
  ylab('Fold Change Lyst / ActB') + 
  #theme_bw(base_size=9) + 
  ggtitle('Lyst fold change by exon')
p1
ggsave(filename='Lyst_fold-change_984-WT_985-KO_qPCR-2022-04-03.png')

print(Dat)


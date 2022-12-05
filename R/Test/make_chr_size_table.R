#############
## 1KG data #
#############
library(data.table)
library(dplyr)

index.1kg <- "/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_index.gz"
index.1kg <- fread(index.1kg)
colnames(index.1kg) <- c("rsid","chr","bp","a1","a2","af1","fp")
head(index.1kg)
dim(index.1kg) #38219238        7

chr_size_1kg <- index.1kg %>% 
  group_by(chr) %>% 
  summarize(p.start=min(bp),q.end=max(bp),num.snp=n())

## find chromosome p arm end bp and q arm start bp
p.end <- rep(NA,22)
q.start <- rep(NA,22)
for(i in 1:22){
  chr <- subset(index.1kg, chr==i)
  bp1 <- chr$bp[1:(length(chr$bp)-1)]
  bp2 <- chr$bp[2:length(chr$bp)]
  if(i %in% c(1:12,16:20)){
    p.end[i] <- bp1[(bp2-bp1)>2000000]
    q.start[i] <- bp2[(bp2-bp1)>2000000]
  } else { # chr 13, 14, 15, 21, 22 only has q arm. (i.e., acrocentric chromosomes https://en.wikipedia.org/wiki/Centromere#Acrocentric)
    p.end[i] <- NA
    q.start[i] <- chr_size_1kg$p.start[i]
    chr_size_1kg$p.start[i] <- NA
  }
  print(i)
}
chr_size_1kg$p.end <- p.end
chr_size_1kg$q.start <- q.start
chr_size_1kg <- chr_size_1kg[,c(1,2,5,6,3,4)]
print(chr_size_1kg, n=22)

chr_size_1kg$num.1mb.p <- round((chr_size_1kg$p.end-chr_size_1kg$p.start)/1000000, digits=2)
chr_size_1kg$num.1mb.q <- round((chr_size_1kg$q.end-chr_size_1kg$q.start)/1000000, digits=2)
head(chr_size_1kg)

write.table(chr_size_1kg, file="/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_chr_size.txt", sep="\t", row.names=FALSE, quote=FALSE)
#read.table(file="/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_chr_size.txt", header=TRUE)

##############
## 33KG data #
##############

index.33kg <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_index.gz"
index.33kg <- fread(index.33kg)
colnames(index.33kg) <- c("rsid","chr","bp","a1","a2","af1","fp")
head(index.33kg)
dim(index.33kg) # 24074297        7


chr_size_33kg <- index.33kg %>% 
  group_by(chr) %>% 
  summarize(p.start=min(bp),q.end=max(bp),num.snp=n())
head(chr_size_33kg)

## find chromosome p arm end bp and q arm start bp
p.end <- rep(NA,22)
q.start <- rep(NA,22)
for(i in 1:22){
  chr <- subset(index.33kg, chr==i)
  bp1 <- chr$bp[1:(length(chr$bp)-1)]
  bp2 <- chr$bp[2:length(chr$bp)]
  if(i %in% c(1:12,16:20)){
    p.end[i] <- bp1[(bp2-bp1)>2000000]
    q.start[i] <- bp2[(bp2-bp1)>2000000]
  } else { # chr 13, 14, 15, 21, 22 only has q arm. (i.e., acrocentric chromosomes https://en.wikipedia.org/wiki/Centromere#Acrocentric)
    p.end[i] <- NA
    q.start[i] <- chr_size_33kg$p.start[i]
    chr_size_33kg$p.start[i] <- NA
  }
  print(i)
}
chr_size_33kg$p.end <- p.end
chr_size_33kg$q.start <- q.start
chr_size_33kg <- chr_size_33kg[,c(1,2,5,6,3,4)]
print(chr_size_33kg, n=22)

chr_size_33kg$num.1mb.p <- round((chr_size_33kg$p.end-chr_size_33kg$p.start)/1000000, digits=2)
chr_size_33kg$num.1mb.q <- round((chr_size_33kg$q.end-chr_size_33kg$q.start)/1000000, digits=2)
print(chr_size_33kg, n=22)

write.table(chr_size_33kg, file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_chr_size.txt", sep="\t", row.names=FALSE, quote=FALSE)



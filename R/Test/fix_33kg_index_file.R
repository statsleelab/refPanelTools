## 33kg_index.gz is not sorted in increasing order of chromosome and base pair poisiton.

#library(bit64)
#.Machine$integer.max

library(data.table)
library(MakeRefPanel)
index.data <- fread("/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_index.gz")
dim(index.data) # 24074297 6
head(index.data)

fpoffset <- fread("/Users/leed13/Desktop/GAUSS/ref/Human/33KG/fpoffset_length.txt")
head(fpoffset)
dim(fpoffset)
sub.dat <- subset(fpoffset, V2!=33216)
head(sub.dat)
dim(sub.dat)

ref.data <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_geno.gz"

geno.norm.info1 <- get_geno_info(as.numeric(index.data$V6[1]), ref.data) #33216
geno.norm.info2 <- get_geno_info(as.numeric(index.data$V6[2]), ref.data) #33216
geno.norm.info3 <- get_geno_info(as.numeric(index.data$V6[3]), ref.data) #33216

nchar(geno.norm.info1)
nchar(geno.norm.info2)
nchar(geno.norm.info3)

af1ref <- fread("/Users/leed13/Desktop/GAUSS/ref/Human/33KG/af1ref.txt")
head(af1ref)
dim(af1ref)

## add af1ref in index.data

index.data$af1ref <- af1ref
head(index.data)
index.data <- index.data[,c(1:5,7,6)]
head(index.data)
index.data <- as.data.frame(index.data)
class(index.data)

## sava as a gz file 
##gzfw = gzfile("/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_index_af1ref.gz","w")
options(scipen = 999) # disable scientific notation
write.table(index.data, gzfw, quote=FALSE, sep=" ", 
            row.names = FALSE, col.names = FALSE)
options(scipen = 0) # restore the default
##close(gzfw)







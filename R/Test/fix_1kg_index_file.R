# Find 1kg data problems and fix.

library(data.table)
library(MakeRefPanel)
index.1kg <- fread("/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_index.gz")
head(index.1kg) 
dim(index.1kg) #38219238 7


af1ref <- fread("/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/af1ref.txt.gz")
head(af1ref)
dim(af1ref)

sum(index.1kg$V6!=af1ref$V1) 
diff <- af1ref$V1 - index.1kg$V6
max(diff) # difference is negligible


fpoffset <- fread("/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/fpoffset_length.txt.gz")
head(fpoffset)
dim(fpoffset)
sum(index.1kg$V7!=fpoffset$V1)


sub.dat <- subset(fpoffset, V2!=1219)
head(sub.dat)
dim(sub.dat)
ref.data <- "/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_geno_af1.gz"

geno.info1 <- get_geno_info(as.numeric(sub.dat$V1[1]), ref.data)
geno.info2 <- get_geno_info(as.numeric(sub.dat$V1[2]), ref.data)
geno.info3 <- get_geno_info(as.numeric(sub.dat$V1[3]), ref.data)
geno.info4 <- get_geno_info(as.numeric(sub.dat$V1[4]), ref.data)
geno.info5 <- get_geno_info(as.numeric(sub.dat$V1[5]), ref.data)
geno.info6 <- get_geno_info(as.numeric(sub.dat$V1[6]), ref.data)
geno.info7 <- get_geno_info(as.numeric(sub.dat$V1[7]), ref.data)
geno.norm.info1 <- get_geno_info(as.numeric(index.1kg$V7[1]), ref.data) #1219
geno.norm.info2 <- get_geno_info(as.numeric(index.1kg$V7[2]), ref.data) #1219
geno.norm.info3 <- get_geno_info(as.numeric(index.1kg$V7[3]), ref.data) #1219

nchar(geno.info1)
nchar(geno.info2)
nchar(geno.info3)
nchar(geno.info4)
nchar(geno.info5)

get_subject_num <- function(geno.info){
  geno.info.vec <- scan(text=geno.info, what="") 
  geno.vec <- geno.info.vec[1:14]
  geno.vec.num <- nchar(geno.vec)
  #af1.vec <- as.numeric(geno.info.vec[15:28])
  geno.vec.num
}
get_subject_num(geno.info1)
get_subject_num(geno.info2)
get_subject_num(geno.info3)
get_subject_num(geno.info4)
get_subject_num(geno.info5)
get_subject_num(geno.info6)
get_subject_num(geno.info7)


## test
library(MakeRefPanel)
fpos <- as.numeric("185606618961528")
ref.data <- "/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_geno_af1.gz"
get_geno_info(fpos, ref.data)

library(data.table)
library(MakeRefPanel)
index.1kg <- fread("/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_index.gz")
head(index.1kg) 
dim(index.1kg) #38219238 7



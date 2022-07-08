#############
## 1KG data #
#############
library(refPanelTools)
ref.data <- "/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_geno_af1.gz"
output <- "/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/fpoffset_length.txt"
indexer(ref.data, output)

# calculate reference allele freq of each SNP in 1KG data
library(refPanelTools)
ref.data <- "/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_geno_af1.gz"
output <- "/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/af1ref.txt" 
num.pops <- 14
cal_af1ref(ref.data, num.pops, output)

##############
## 33KG data #
##############
library(refPanelTools)
ref.data <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_geno.gz"
output <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/fpoffset_length.txt" 
indexer(ref.data, output)


# calculate reference allele freq of each SNP in 33KG data
library(refPanelTools)
ref.data <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_geno.gz"
output <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/af1ref.txt" 
num.pops <- 29
cal_af1ref(ref.data, num.pops, output)


#####################################################
## Make a reference panel file for each chromosome ##
#####################################################
library(data.table)
library(refPanelTools)
#chr.num <- 20
num.pops <- 29
ref.index.file <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_index.gz"
ref.data.file <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_geno.gz"

ref.index <- fread(ref.index.file)
head(ref.index)

for(chr.num in 1:9){
  chr.index <- subset(ref.index, V2==chr.num)
  head(chr.index)
  dim(chr.index)
  
  # extract chr genotype data
  data.output <- paste0("/Users/leed13/Desktop/GAUSS/ref/Human/33KG/chr_data/33kg_chr",chr.num,"_geno")
  extract_chr_data(chr.num, num.pops, 
                   ref.index.file, ref.data.file, 
                   data.output)
  system(paste0("bgzip ",data.output))
  
  # do indexing genotype data
  chr.data.file <- paste0(data.output,".gz")
  fp.output <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/chr_data/fpoffset_length_chr.txt" 
  indexer(chr.data.file, fp.output)
  
  # read fp offsets  
  chr.fp <- fread(fp.output)
  head(chr.index)
  head(chr.fp)
  chr.index$V7 <- chr.fp$V1
  head(chr.index)
  
  chr.index.file <- paste0("/Users/leed13/Desktop/GAUSS/ref/Human/33KG/chr_data/33kg_chr",chr.num,"_index")
  fwrite(chr.index, file=chr.index.file, quote=FALSE, sep=" ",row.names = FALSE, col.names = FALSE)
  system(paste0("bgzip ",chr.index.file))
  system(paste0("rm ",fp.output))
}

#############################################################
## Extract genotype data of a user-specifed genomic region ##
#############################################################
library(refPanelTools)
chr.num <- 14
start.bp <- 104000000
end.bp   <- 104200000
num.pops <- 29
ref.index.file <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_index.gz"
ref.data.file <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_geno.gz"

# extract genotype data of a user-specified genomic region 
data.output <- paste0("/Users/leed13/Desktop/GAUSS/ref/Human/33KG/chr_data/33kg_chr",chr.num,"_reg_geno")
extract_reg_data(chr.num, start.bp, end.bp, num.pops, 
                 ref.index.file, ref.data.file, 
                 data.output)
#system(paste0("bgzip ",data.output))










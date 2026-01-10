#############
## 1KG data #
#############
library(refPanelTools)

fpath <- "/Users/leed13/Documents/My_Project/GAUSS_Reference_Panels/Human"

ref.data <- paste0(fpath,"/1KG/v0.1.0/1kg_geno_af1.gz")
output   <- paste0(fpath,"/1KG/v0.1.0/fpoffset_length.txt")
indexer(ref.data, output)

# calculate reference allele freq of each SNP in 1KG data
library(refPanelTools)
ref.data <- paste0(fpath,"/1KG/v0.1.0/1kg_geno_af1.gz")
output   <- paste0(fpath,"/1KG/v0.1.0/af1ref.txt") 
num.pops <- 14
cal_af1ref(ref.data, num.pops, output)

##############
## 33KG data #
##############
library(refPanelTools)
ref.data <- paste0(fpath,"/33KG/33kg_geno.gz")
output   <- paste0(fpath,"/33KG/fpoffset_length.txt") 
indexer(ref.data, output)


# calculate reference allele freq of each SNP in 33KG data
library(refPanelTools)
ref.data <- paste0(fpath,"/33KG/33kg_geno.gz")
output   <- paste0(fpath,"/33KG/af1ref.txt") 
num.pops <- 29
cal_af1ref(ref.data, num.pops, output)


#####################################################
## Make a reference panel file for each chromosome ##
#####################################################
library(data.table)
library(refPanelTools)
#chr.num <- 20
num.pops <- 29
ref.index.file <- paste0(fpath,"/33KG/33kg_index.gz")
ref.data.file  <- paste0(fpath,"/33KG/33kg_geno.gz")

ref.index <- fread(ref.index.file)
head(ref.index)

for(chr.num in 1:9){
  chr.index <- subset(ref.index, V2==chr.num)
  head(chr.index)
  dim(chr.index)
  
  # extract chr genotype data
  data.output <- paste0(fpath,"/33KG/chr_data/33kg_chr",chr.num,"_geno")
  extract_chr_data(chr.num, num.pops, 
                   ref.index.file, ref.data.file, 
                   data.output)
  system(paste0("bgzip ",data.output))
  
  # do indexing genotype data
  chr.data.file <- paste0(data.output,".gz")
  fp.output <- paste0(fpath,"/33KG/chr_data/fpoffset_length_chr.txt") 
  indexer(chr.data.file, fp.output)
  
  # read fp offsets  
  chr.fp <- fread(fp.output)
  head(chr.index)
  head(chr.fp)
  chr.index$V7 <- chr.fp$V1
  head(chr.index)
  
  chr.index.file <- paste0(fpath,"/33KG/chr_data/33kg_chr",chr.num,"_index")
  fwrite(chr.index, file=chr.index.file, quote=FALSE, sep=" ",row.names = FALSE, col.names = FALSE)
  system(paste0("bgzip ",chr.index.file))
  system(paste0("rm ",fp.output))
}

#############################################################
## Extract genotype data of a user-specified genomic region ##
#############################################################
library(refPanelTools)
chr.num <- 14
start.bp <- 104000000
end.bp   <- 104200000
num.pops <- 29
ref.index.file <- paste0(fpath,"/33KG/33kg_index.gz")
ref.data.file  <- paste0(fpath,"/33KG/33kg_geno.gz")

# extract genotype data of a user-specified genomic region 
data.output <- paste0(fpath,"/33KG/chr_data/33kg_chr",chr.num,"_reg_geno")
extract_reg_data(chr.num, start.bp, end.bp, num.pops, 
                 ref.index.file, ref.data.file, 
                 data.output)
#system(paste0("bgzip ",data.output))


########################################################################
## Extract genotype and allele frequency data of specific populations ##
########################################################################
library(refPanelTools)
chr.num <- 22
ref.index.file <- paste0(fpath,"/33KG/33kg_index.gz")
ref.data.file  <- paste0(fpath,"/33KG/33kg_geno.gz")
ref.desc.file  <- paste0(fpath,"/33KG/33kg_pop_desc.txt")
pop.vec <- c("CCE","CEU")

# extract genotype data of a user-specified genomic region 
data.output <- paste0(fpath,"/33KG/test/33kg_chr",chr.num,"_CCE_CEU.txt")
extract_chr_pop_data(chr.num, 
                     pop.vec, 
                     ref.index.file, 
                     ref.data.file,
                     ref.desc.file,
                     data.output)


library(data.table)
chr22 <- fread(data.output, sep=" ", colClasses = 'character')
dim(chr22)


####################################################################
## Extract allele frequency data of all ethnic groups populations ##
####################################################################
library(refPanelTools)
chr.num <- 22
ref.index.file <- paste0(fpath,"/33KG/33kg_index.gz")
ref.data.file  <- paste0(fpath,"/33KG/33kg_geno.gz")
ref.desc.file  <- paste0(fpath,"/33KG/33kg_pop_desc.txt")

# extract af1 of all ethnic groups 
data.output <- paste0(fpath,"/33KG/test/33kg_chr",chr.num,"_af1.txt")
extract_all_af1(chr.num, 
                ref.index.file, 
                ref.data.file,
                ref.desc.file,
                data.output)


library(data.table)
chr22 <- fread(data.output, sep=" ", colClasses = 'character')
dim(chr22)


##########################################################
## Simulate SNP allele frequencies of multi-ethnic GWAS ##
##########################################################
library(refPanelTools)
chr.num <- 22
ref.index.file <- paste0(fpath,"/33KG/33kg_index.gz")
ref.data.file <-  paste0(fpath,"/33KG/33kg_geno.gz")
ref.desc.file <-  paste0(fpath,"/33KG/33kg_pop_desc.txt")

pop.vec <- c("CCE","CEU")
pop.num.vec <- c(5000,5000)
#pop.num.vec <- c(90000,10000)

#pop.vec <- c("BEB","CLM")
#pop.num.vec <- c(10,10)

data.output <- paste0(fpath,"/33KG/test/33kg_chr",chr.num,"_",pop.vec[1],"_",pop.vec[2],"_",pop.num.vec[1],"_",pop.num.vec[2],"_sim_af1_sim_z.txt")
data.output
simulate_af1_z(chr.num, 
             pop.vec, 
             pop.num.vec,
             ref.index.file, 
             ref.data.file,
             ref.desc.file,
             data.output)

library(data.table)
chr22 <- fread(data.output, sep=" ", header=TRUE)
head(chr22)
dim(chr22)


hist(chr22$sim_z)
plot(chr22$sim_af1, chr22$sim_z)
#plot(chr22$sim_af1, chr22$beta1)
#plot(chr22$sim_af1, chr22$beta0)
#plot(chr22$sim_af1, chr22$std_err)
#plot(chr22$sim_af1, chr22$mse)
#plot(chr22$sim_af1, chr22$sxy)
#plot(chr22$sim_af1, chr22$sxx)




library(gauss)

wgt.df1 <- afmix(input_file=data.output,
                 reference_index_file = ref.index.file,
                 reference_data_file = ref.data.file,
                 reference_pop_desc_file = ref.desc.file,
                 interval=10)

wgt.df2 <- cpw2(input_file=data.output,
                reference_index_file = ref.index.file,
                reference_data_file = ref.data.file,
                reference_pop_desc_file = ref.desc.file,
                interval=10)


wgt.df1
wgt.df2

new.wgt1 <- subset(wgt.df1, wgt>0.01)
new.wgt2 <- subset(wgt.df2, wgt>0.01)
new.wgt1
new.wgt2
sum(new.wgt1$wgt)
sum(new.wgt2$wgt)


#pop   wgt
#1 CCE 0.407
#2 CCS 0.054
#3 CEU 0.391
#4 CSE 0.044
#8 ORK 0.105

library(tidyverse)
wgt.df1 %>%
  group_by(sup.pop) %>%
  summarise(wgt=sum(wgt), .groups="drop")

# A tibble: 5 × 2
#sup.pop   wgt
#<chr>   <dbl>
#1 AFR     0.002
#2 AMR     0.003
#3 ASN     0.505
#4 EUR     0.508
#5 SAS     0.002



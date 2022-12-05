#############################################
## Lift over using MungeSumstats R package ##
#############################################

## 1. Install the most recent R version > 4.2 from https://www.r-project.org/

## 2. Install Bioconductor package 3.15
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

## 3. Install R package MungeSumstats
## https://bioconductor.org/packages/release/bioc/html/MungeSumstats.html
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MungeSumstats")

## 4. Install GRCh38 reference data
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
## 5. Install GRCh37 reference data
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

## 5. Load MungeSumstats
library(MungeSumstats)

## 6. Read GWAS summary stats you want to lift over

## Example1 convert hg19 to hg38. 
if(FALSE){
sumstats_dt <- MungeSumstats::formatted_example()
sumstats_dt_hg38 <- MungeSumstats::liftover(sumstats_dt = sumstats_dt, 
                                            ref_genome = "hg19",
                                            convert_ref_genome = "hg38")
}

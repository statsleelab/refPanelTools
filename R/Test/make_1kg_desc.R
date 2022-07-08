pop.1kg.dat <- data.frame(
  Population_Abbreviation=factor(c("ASW","CEU","CHB","CHS","CLM","FIN","GBR","IBS","JPT","LWK","MXL","PUR","TSI","YRI")),
  Number_of_Subject=c(61,85,97,100,60,93,89,14,89,97,66,55,98,88),
  Super_Population=factor(c("AFR","EUR","ASN","ASN","AMR","EUR","EUR","EUR","ASN","AFR","AMR","AMR","EUR","AFR")),
  Population_Description=c("African Ancestry in Southwest US",
                           "Utah residents with Northern and Western European ancestry",
                           "Han Chinese in Beijing, China",
                           "Southern Han Chinese",
                           "Colombian in Medellin, Colombia",
                           "Finnish in Finland",
                           "British in England and Scotlant",
                           "Iberian populations in Spain",
                           "Japanese in Tokyo, Japan",
                           "Luhya in Wenbuye, Kenya",
                           "Mexican Ancestry from Los Angeles, USA",
                           "Puerto Rican in Puerto Rico",
                           "Toscani in Italia",
                           "Yoruba in Ibadan, Nigeria")
)
pop.1kg.dat

write.table(pop.1kg.dat, file="/Users/leed13/Desktop/GAUSS/ref/Human/1KG/v0.1.0/1kg_pop_list.txt", sep="\t", row.names=FALSE, quote=FALSE)

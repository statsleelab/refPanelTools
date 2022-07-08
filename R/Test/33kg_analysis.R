data <- read.table("/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_geno_10lines.txt", sep="", colClasses="character", header=FALSE, blank.lines.skip=FALSE)
dim(data)
geno <- data[,1:29]
af1 <- data[,30:58]
nchar_vec <- nchar(geno[1,])
nchar_vec

pop.33kg.dat <- data.frame(
  Population_Abbreviation=factor(c("ACB","ASW","BEB","CCE","CCS","CDX","CEU","CLM","CNE","CSE","ESN","FIN","GBR","GIH","GWD","IBS"
                            ,"ITU","JPT","KHV","LWK","MSL","MXL","ORK","PEL","PJL","PUR","STU","TSI","YRI")),
  Number_of_Subject=c(164,162,86,3409,2613,95,6360,98,2330,2020,140,3529,2020,110,113,1309,95,107,226,99,87,187,5772,110,121,
                      138,110,1291,52),
  Super_Population=factor(c("AFR","AFR","SAS","ASN","ASN","ASN","EUR","AMR","ASN","ASN","AFR","EUR","EUR","SAS","AFR","EUR","SAS",
                     "ASN","ASN","AFR","AFR","AMR","EUR","AMR","SAS","AMR","SAS","EUR","AFR")),
  Population_Description=c("African Caribbeans in Barbados","African Ancestry in Southwest US","Bengali from Bangladesh",
                           "China Central East","China Central South","Chinese Dai in Xishuangbanna, China",
                           "Utah residents with Northern and Western European ancestry","Colombians from Medellin, Colombia",
                           "China North East","China South-East","Esan in Nigeria","Finnish in Finland","British in England and Scotland",
                           "Gujarati Indian from Houston, Texas","Gambian in Western Divisions in the Gambia","Iberian Population in Spain",
                           "Indian Telugu from the UK","Japanese in Tokyo, Japan","Kinh in Ho Chi Minh City, Vietnam","Luhya in Webuye, Kenya",
                           "Mende in Sierra Leone","Mexican Ancestry from Los Angeles, USA","Orkney Island study","Peruvians from Lima, Peru",
                           "Punjabi from Lahore, Pakistan","Puerto Rican in Puerto Rico","Sri Lankan Tamil from the UK","Toscani in Italia",
                           "Yoruba in Ibadan, Nigeria")
)
pop.33kg.dat
pop.33kg.dat$Number_of_Subject==nchar_vec

write.table(pop.33kg.dat, file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg.pop.dat.txt", sep="\t", row.names=FALSE, quote=FALSE)







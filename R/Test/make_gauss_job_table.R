library(data.table)
file <- "/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_chr_size.txt"
chr.table <- fread(file, header=TRUE)

new.table <- data.frame(chr=integer(), arm=character(), start.bp=integer(), end.bp=integer(), num.1mb=numeric())

for(i in 1:nrow(chr.table)){
  ## p arm
  if(!is.na(chr.table$p.start[i])){
    chr <- chr.table$chr[i]
    start.bp <- chr.table$p.start[i]
    end.bp <- chr.table$p.end[i]
    num.1mb <- chr.table$num.1mb.p[i]
    new.table <- rbind(new.table, data.frame(chr=chr, arm="p", start.bp=start.bp, end.bp=end.bp, num.1mb=num.1mb))
  }
  ## q arm
  if(!is.na(chr.table$q.start[i])){
    chr <- chr.table$chr[i]
    start.bp <- chr.table$q.start[i]
    end.bp <- chr.table$q.end[i]
    num.1mb <- chr.table$num.1mb.q[i]
    new.table <- rbind(new.table, data.frame(chr=chr, arm="q", start.bp=start.bp, end.bp=end.bp, num.1mb=num.1mb))
  }
}

write.table(new.table, file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_chr_size_long.txt", sep="\t", row.names=FALSE, quote=FALSE)

##########################
# Make cluster job table #
##########################

new.table <- read.table(file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_chr_size_long.txt", sep="\t", header = TRUE)

win.size <- 100000
#win.size <- 250000
#win.size <- 1000000
job.table <- data.frame(chr=integer(), arm=character(), start.bp=integer(), end.bp=integer())
for(i in 1:nrow(new.table)){
  chr <- new.table$chr[i]
  arm <- new.table$arm[i]
  start.bp <- new.table$start.bp[i]
  end.bp <- new.table$end.bp[i]
  j <- 1
  while(TRUE){
    win.sbp <- start.bp + (win.size)*(j-1)
    win.ebp <- start.bp + (win.size*j) - 1
    if(win.ebp > end.bp){
      job.table <- rbind(job.table, data.frame(chr=chr, arm=arm, start.bp=win.sbp, end.bp=end.bp))
      break;
    }
    job.table <- rbind(job.table, data.frame(chr=chr, arm=arm, start.bp=win.sbp, end.bp=win.ebp))
    j <- j+1
  }
}
job.table$win.size <- job.table$end.bp-job.table$start.bp + 1
job.table <- format(job.table, scientific=FALSE)
write.table(job.table, file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_100kb_job_table.txt", sep="\t", row.names=FALSE, quote=FALSE)# 27043 5
#write.table(job.table, file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_250kb_job_table.txt", sep="\t", row.names=FALSE, quote=FALSE) # 10829 5
#write.table(job.table, file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_1mb_job_table.txt", sep="\t", row.names=FALSE, quote=FALSE) # 2724 5


t1mb <- read.table(file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_1mb_job_table.txt", sep="\t", header=TRUE)
t250kb <- read.table(file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_250kb_job_table.txt", sep="\t", header=TRUE)
t100kb <- read.table(file="/Users/leed13/Desktop/GAUSS/ref/Human/33KG/33kg_100kb_job_table.txt", sep="\t", header=TRUE)

